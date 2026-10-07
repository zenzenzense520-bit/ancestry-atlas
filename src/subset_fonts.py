"""Rebuild renamed CJK font subsets for all report/source text."""
import json
from pathlib import Path
from fontTools.ttLib import TTFont
from fontTools import subset

ROOT = Path(__file__).resolve().parents[1]

def main():
    package = ROOT/'references/fontsource/package'
    text = ''.join(p.read_text() for p in (ROOT/'src').glob('*') if p.is_file() and p.suffix in ['.html','.py','.js','.css'])
    text += ''.join(p.read_text() for p in (ROOT/'results').glob('*.json'))
    wanted = {ord(c) for c in text if ord(c)>127}
    outstanding = wanted.copy(); css=[]
    target = ROOT/'src/fonts';target.mkdir(exist_ok=True)
    for path in sorted((package/'files').glob('*.woff2')):
        font = TTFont(path)
        chars = set(font.getBestCmap()) & outstanding
        if not chars:continue
        options = subset.Options();options.flavor='woff2'
        sub = subset.Subsetter(options=options);sub.populate(unicodes=chars);sub.subset(font)
        for record in font['name'].names:
            if record.nameID in [1,4,6,16]:
                value = 'AtlasSerifSC' if record.nameID==6 else 'Atlas Serif SC'
                record.string = value.encode(record.getEncoding())
        font.flavor='woff2';font.save(target/path.name)
        css.append("@font-face{font-family:'Atlas Serif SC';font-style:normal;font-weight:200 900;font-display:swap;src:url('fonts/"+path.name+"') format('woff2');unicode-range:"+','.join(f'U+{cp:X}' for cp in sorted(chars))+';}')
        outstanding -= chars
        if not outstanding:break
    # Superscript mathematical glyphs may fall back to system fonts. CJK must be covered.
    cjk_missing = [cp for cp in outstanding if 0x3400<=cp<=0x9fff]
    assert not cjk_missing, cjk_missing
    (ROOT/'src/fonts.css').write_text('\n'.join(css))
    (target/'OFL.txt').write_text((package/'LICENSE').read_text())
    (ROOT/'results/font_coverage.json').write_text(json.dumps({'required_non_ascii':len(wanted),'covered':len(wanted-outstanding),
        'uncovered_codepoints':[f'U+{c:X}' for c in sorted(outstanding)],'cjk_missing':cjk_missing,'font_family':'Atlas Serif SC'}))
    print(f'Subset {len(wanted-outstanding)} characters in {len(css)} fonts')

if __name__=='__main__':main()
