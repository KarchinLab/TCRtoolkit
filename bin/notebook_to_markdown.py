#!/usr/bin/env python3
"""Turn an executed Quarto notebook (.quarto_ipynb, kept with keep-ipynb) into plain markdown
that another Quarto document can {{< include >}}, whatever its engine.

The Discovery Brief (Python) is shown inside the Master Summary (R/knitr). One Quarto document
runs one engine, and Quarto's own {{< embed >}} needs a newer Quarto than the single-cell
container ships, so the executed notebook is flattened here instead:

  * markdown and raw-HTML cells are copied as they are (the YAML front matter is dropped);
  * code is never shown; a cell's outputs are copied unless it is `include: false` or
    `output: false`. `output: asis` text is copied verbatim, other printed text goes in a
    code block, HTML in a raw HTML block and images as inline data URIs;
  * Plotly figures load plotly.js through RequireJS from a CDN. The host has no RequireJS,
    so plotly.js is inlined once at the top, a minimal `require` hands figures the global
    Plotly, and the notebook's own RequireJS loader cell is dropped.

Usage: notebook_to_markdown.py <in.quarto_ipynb> <out.md>
"""
import json
import re
import sys


def cell_options(source):
    opts = {}
    for line in source.splitlines():
        m = re.match(r'\s*#\|\s*([\w-]+)\s*:\s*(.*?)\s*$', line)
        if m:
            opts[m.group(1)] = m.group(2).strip('"\'')
        elif line.strip() and not line.strip().startswith('#'):
            break
    return opts


def raw_html(html):
    # A fence longer than any backtick run inside keeps the block intact.
    fence = '`' * max(3, max((len(r) for r in re.findall(r'`+', html)), default=0) + 1)
    return f"\n{fence}{{=html}}\n{html}\n{fence}\n"


def is_plotly_loader(html):
    return 'window.PlotlyConfig' in html and 'plotly-graph-div' not in html


def output_markdown(out, asis):
    kind = out.get('output_type')
    if kind == 'stream':
        if out.get('name') == 'stderr':
            return ''
        text = ''.join(out.get('text', ''))
        return text if asis else (f"\n```\n{text.rstrip()}\n```\n" if text.strip() else '')
    if kind == 'error':
        return f"\n```\n{out.get('ename', '')}: {out.get('evalue', '')}\n```\n"
    data = out.get('data', {})
    if 'text/markdown' in data:
        return '\n' + ''.join(data['text/markdown']) + '\n'
    if 'text/html' in data:
        html = ''.join(data['text/html'])
        return '' if is_plotly_loader(html) else raw_html(html)
    for mime in ('image/png', 'image/jpeg'):
        if mime in data:
            b64 = ''.join(data[mime]).replace('\n', '')
            return raw_html(f'<img src="data:{mime};base64,{b64}" style="max-width:100%;">')
    if 'image/svg+xml' in data:
        return raw_html(''.join(data['image/svg+xml']))
    if 'text/plain' in data and not asis:
        return f"\n```\n{''.join(data['text/plain']).rstrip()}\n```\n"
    return ''


def plotly_preamble():
    try:
        from plotly.offline import get_plotlyjs
    except ImportError:
        return ''
    shim = ("window.require = window.require || function (deps, cb) { cb(window.Plotly); };\n"
            "window.require.undef = window.require.undef || function () {};")
    return raw_html(f'<script type="text/javascript">{get_plotlyjs()}</script>\n'
                    f'<script type="text/javascript">{shim}</script>')


def main(src, dst):
    nb = json.load(open(src))
    parts = []
    has_plotly = False
    for i, cell in enumerate(nb['cells']):
        source = ''.join(cell['source']) if isinstance(cell['source'], list) else cell['source']
        ctype = cell['cell_type']
        if ctype == 'markdown':
            if i == 0:
                source = re.sub(r'\A---\n.*?\n---\n', '', source, flags=re.S)
            parts.append('\n' + source.strip('\n') + '\n')
        elif ctype == 'raw':
            fmt = cell.get('metadata', {}).get('raw_mimetype', '')
            parts.append(raw_html(source) if 'html' in fmt else '\n' + source + '\n')
        elif ctype == 'code':
            opts = cell_options(source)
            if opts.get('include') == 'false' or opts.get('output') == 'false':
                continue
            asis = opts.get('output') == 'asis'
            for out in cell.get('outputs', []):
                md = output_markdown(out, asis)
                has_plotly = has_plotly or 'plotly-graph-div' in md
                parts.append(md)
    body = '\n'.join(p for p in parts if p)
    with open(dst, 'w') as fh:
        if has_plotly:
            fh.write(plotly_preamble())
        fh.write(body + '\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
