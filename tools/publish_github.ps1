param(
    [string]$Repository = 'zenzenzense520-bit/ancestry-atlas'
)
$ErrorActionPreference = 'Stop'
$ProjectRoot = if ($PSScriptRoot) { Split-Path -Parent $PSScriptRoot } elseif ($env:ATLAS_PUBLIC_ROOT) { $env:ATLAS_PUBLIC_ROOT } else { throw 'Cannot locate the source project.' }
$ExpectedOwner = 'zenzenzense520-bit'
if ($Repository -ne "$ExpectedOwner/ancestry-atlas") {
    throw 'This publisher is scoped to zenzenzense520-bit/ancestry-atlas.'
}
foreach ($command in @('git', 'gh')) {
    if (-not (Get-Command $command -ErrorAction SilentlyContinue)) {
        throw "Install $command from its official website, then run this publisher again."
    }
}
function Read-Native {
    param([string]$Command, [string[]]$NativeArguments)
    $previous = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        $output = @(& $Command @NativeArguments 2>$null)
        $code = $LASTEXITCODE
        return [pscustomobject]@{ Output = ($output -join "`n"); ExitCode = $code }
    } finally { $ErrorActionPreference = $previous }
}
$auth = Read-Native -Command 'gh' -NativeArguments @('auth', 'status', '--hostname', 'github.com')
if ($auth.ExitCode -ne 0) {
    throw 'GitHub CLI is not signed in. Run: gh auth login --hostname github.com --git-protocol https --web'
}
$login = Read-Native -Command 'gh' -NativeArguments @('api', 'user', '--jq', '.login')
if ($login.ExitCode -ne 0 -or $login.Output.Trim() -ne $ExpectedOwner) {
    throw "GitHub CLI must be signed in as $ExpectedOwner."
}
function Invoke-Git {
    param([string[]]$GitArguments)
    & git @GitArguments
    if ($LASTEXITCODE -ne 0) {
        throw "Git failed: $($GitArguments -join ' ')"
    }
}
Push-Location $ProjectRoot
try {
    if (-not (Test-Path -LiteralPath (Join-Path $ProjectRoot '.git'))) {
        Invoke-Git -GitArguments @('init', '--initial-branch=main')
    }
    $branch = & git branch --show-current
    if ($LASTEXITCODE -ne 0 -or $branch.Trim() -ne 'main') {
        throw 'This publisher requires the main branch.'
    }
    $name = & git config user.name
    if (-not $name) { Invoke-Git -GitArguments @('config', 'user.name', $ExpectedOwner) }
    $email = & git config user.email
    if (-not $email) {
        Invoke-Git -GitArguments @('config', 'user.email', '239009337+zenzenzense520-bit@users.noreply.github.com')
    }
    Invoke-Git -GitArguments @('add', '--', 'README.md', 'LICENSE', 'requirements.txt', '.gitignore', 'src', 'tests', 'docs', 'tools', '.github', 'Publish-GitHub.cmd')
    $tracked = & git ls-files
    if ($LASTEXITCODE -ne 0) { throw 'Cannot inspect publication scope.' }
    $denied = $tracked | Where-Object {
        $_ -match '(^|/)(input|upload|results|report|figures|qa|references|runtime|deliverables|analysis-runtime|research|\.venv|venv|__pycache__)(/|$)' -or
        $_ -match '\.(npz|npy|pgen|pvar|psam|bed|bim|fam|vcf|bcf|raw|zip|zst|pyc)$' -or
        $_ -match '(^|/)\.env($|\.)'
    }
    if ($denied) { throw "Private files in publication scope: $($denied -join ', ')" }
    & git diff --cached --quiet
    $changes = $LASTEXITCODE
    if ($changes -eq 1) {
        Invoke-Git -GitArguments @('commit', '-m', 'Publish Ancestry Atlas v0.5 source and Windows setup')
    } elseif ($changes -ne 0) {
        throw 'Cannot inspect staged changes.'
    }
    & gh auth setup-git --hostname github.com
    if ($LASTEXITCODE -ne 0) { throw 'GitHub CLI could not configure Git authentication.' }
    $view = Read-Native -Command 'gh' -NativeArguments @('repo', 'view', $Repository, '--json', 'visibility', '--jq', '.visibility')
    $exists = $view.ExitCode -eq 0
    $originResult = Read-Native -Command 'git' -NativeArguments @('remote', 'get-url', 'origin')
    $origin = $originResult.Output
    if ($origin -and $origin.Trim() -notin @("https://github.com/$Repository.git", "https://github.com/$Repository", "git@github.com:$Repository.git")) {
        throw 'The existing origin points to a different repository.'
    }
    if ($exists) {
        if ($view.Output.Trim() -ne 'PUBLIC') { throw 'The existing repository is not public; its visibility was not changed.' }
        if (-not $origin) { Invoke-Git -GitArguments @('remote', 'add', 'origin', "https://github.com/$Repository.git") }
        Invoke-Git -GitArguments @('push', '--set-upstream', 'origin', 'main')
    } else {
        $createArguments = @('repo', 'create', $Repository, '--public', '--source', $ProjectRoot, '--push', '--description', 'Auditable local genotype QC and reference PCA with offline ancestry reports')
        if (-not $origin) { $createArguments += @('--remote', 'origin') }
        & gh @createArguments
        if ($LASTEXITCODE -ne 0) { throw 'Repository creation or upload failed. No force push was attempted.' }
    }
    $remoteState = & gh repo view $Repository --json visibility,url
    if ($LASTEXITCODE -ne 0) { throw 'Upload finished but repository verification failed.' }
    $verified = $remoteState | ConvertFrom-Json
    if ($verified.visibility -ne 'PUBLIC') { throw 'The resulting repository is not public.' }
    & git ls-remote --exit-code origin refs/heads/main
    if ($LASTEXITCODE -ne 0) { throw 'Cannot verify the remote main branch.' }
    Write-Host "Public repository: $($verified.url)"
    Write-Host 'Only the source project was uploaded. The sibling research folder was not uploaded.'
} finally {
    Pop-Location
}
