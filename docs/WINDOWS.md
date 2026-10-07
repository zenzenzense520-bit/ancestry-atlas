# Windows / D 盘

D 盘迁移包解压后，双击根目录的 `Install-D-Drive.cmd`。它会创建 `D:\AncestryAtlas`。同名文件夹已经存在时，使用编号新文件夹，不覆盖已有项目，也不删除下载文件。

- `ancestry-atlas\`：公开源码、报告模板、许可、测试和发布脚本。
- `research\`：本地研究工作区，含原始输入、规范化调用、已保存模型、QC、参考调用交集、报告和恢复记录。
- `research\report\Sample01_reference_report_v05.html`：可离线打开的个人研究报告。
- 大型完整公共参考需要按 README 重新下载；迁移包包含既有分析的参考交集。

公开上传在 `ancestry-atlas\` 内双击 `Publish-GitHub.cmd`。需要官方 [Git](https://git-scm.com/downloads/win) 与 [GitHub CLI](https://cli.github.com/)；GitHub CLI 必须登录 `zenzenzense520-bit`。如未登录，在终端运行：

```powershell
gh auth login --hostname github.com --git-protocol https --web
```

发布脚本新建公开的 `zenzenzense520-bit/ancestry-atlas` 并上传源码。它核对账号、提交范围、公开状态与远程分支，使用普通推送。同名仓库已经存在时，只接受公开仓库；遇到不同远程地址或历史冲突时停止。

可在个人研究目录中复核保存的模型：

```powershell
cd D:\AncestryAtlas\research
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy==2.3.5
$env:OPENBLAS_NUM_THREADS = '2'
$env:OMP_NUM_THREADS = '2'
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_*.py' -v
```

完整重跑的数据下载、绘图和解码还需安装 `requirements.txt`。分析代码此前在 Linux / Python 3.12 验证；Windows 脚本已静态复核，未在用户电脑上执行。

安装包中记录 SHA-256。个人研究文件仅用于本地分析；公开发布脚本只在独立源码目录工作。
