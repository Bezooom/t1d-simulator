# -*- coding: utf-8 -*-
"""
Publication-Quality Manuscript Formatter & PDF / PNG Generator
Converts Markdown manuscript (bioRxiv format) to styled HTML, renders math via KaTeX,
generates PDF via Playwright, and extracts page screenshots for visual verification.
"""

import os
import sys
import re
import markdown
from playwright.sync_api import sync_playwright
import fitz  # PyMuPDF

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS_DIR = os.path.join(BASE_DIR, "docs")
REPORTS_DIR = os.path.join(BASE_DIR, "reports", "manuscript_screenshots")

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <title>bioRxiv Manuscript - T1D Digital Twin</title>
    <!-- KaTeX CSS & JS for LaTeX Math rendering -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.css">
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.js"></script>
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/contrib/auto-render.min.js" 
            onload="renderMathInElement(document.body, {{delimiters: [
                {{left: '$$', right: '$$', display: true}},
                {{left: '$', right: '$', display: false}}
            ]}});"></script>
    <style>
        @page {{
            size: A4;
            margin: 15mm;
            @bottom-right {{
                content: "Page " counter(page) " of " counter(pages);
                font-family: Arial, sans-serif;
                font-size: 9pt;
                color: #666;
            }}
        }}

        body {{
            font-family: 'Georgia', 'Times New Roman', serif;
            font-size: 10.5pt;
            line-height: 1.5;
            color: #1a1a1a;
            margin: 0 auto;
            padding: 10px 20px;
            background-color: #ffffff;
        }}

        /* Header / Banner */
        .biorxiv-header {{
            border-bottom: 3px solid #8b0000;
            padding-bottom: 6px;
            margin-bottom: 18px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .biorxiv-logo {{
            font-family: 'Helvetica Neue', Arial, sans-serif;
            font-size: 22pt;
            font-weight: bold;
            color: #8b0000;
            letter-spacing: -0.5px;
        }}
        .biorxiv-sub {{
            font-family: Arial, sans-serif;
            font-size: 8.5pt;
            color: #555;
            text-transform: uppercase;
            letter-spacing: 1px;
            text-align: right;
        }}

        /* Title & Authors */
        h1 {{
            font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
            font-size: 17pt;
            font-weight: 700;
            line-height: 1.28;
            color: #111;
            margin-top: 5px;
            margin-bottom: 12px;
        }}

        .meta-info {{
            font-family: Arial, sans-serif;
            font-size: 9pt;
            color: #333;
            background: #f8f9fa;
            border-left: 4px solid #8b0000;
            padding: 8px 12px;
            margin-bottom: 18px;
            border-radius: 3px;
        }}

        /* Abstract Box - Prevent Awkward Page Break */
        .abstract-box {{
            font-family: Arial, sans-serif;
            font-size: 9.5pt;
            line-height: 1.55;
            background-color: #f0f4f8;
            border: 1px solid #d0d7de;
            border-radius: 6px;
            padding: 14px 18px;
            margin-bottom: 22px;
            page-break-inside: avoid;
            break-inside: avoid;
        }}
        .abstract-box h2 {{
            font-size: 11pt;
            text-transform: uppercase;
            color: #0969da;
            margin-top: 0;
            margin-bottom: 8px;
            letter-spacing: 0.5px;
            border-bottom: None;
            padding-bottom: 0;
        }}

        /* Headings */
        h2 {{
            font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
            font-size: 12.5pt;
            font-weight: bold;
            color: #0969da;
            border-bottom: 1px solid #eaecef;
            padding-bottom: 4px;
            margin-top: 20px;
            margin-bottom: 10px;
            page-break-after: avoid;
        }}

        h3 {{
            font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
            font-size: 10.5pt;
            font-weight: bold;
            color: #24292f;
            margin-top: 15px;
            margin-bottom: 6px;
            page-break-after: avoid;
        }}

        /* Scientific Figures & Captions */
        .figure-container {{
            margin: 18px auto;
            text-align: center;
            page-break-inside: avoid;
            break-inside: avoid;
        }}

        .figure-container img {{
            max-width: 60%;
            height: auto;
            display: block;
            margin: 0 auto;
            border: 1px solid #d1d5da;
            border-radius: 4px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }}

        .graphical-abstract-img {{
            max-width: 70% !important;
        }}

        .figure-caption {{
            font-family: Arial, sans-serif;
            font-size: 8.5pt;
            color: #444;
            max-width: 82%;
            margin: 8px auto 0 auto;
            line-height: 1.35;
            text-align: justify;
        }}

        /* Formulas & Code */
        .katex-display {{
            margin: 10px 0 !important;
            padding: 6px;
            background: #fafafa;
            border-radius: 4px;
            page-break-inside: avoid;
        }}

        code {{
            font-family: 'Courier New', Courier, monospace;
            background: #f6f8fa;
            padding: 2px 4px;
            border-radius: 3px;
            font-size: 9pt;
        }}

        /* References */
        ol {{
            font-size: 9pt;
            color: #333;
            padding-left: 18px;
        }}
        li {{
            margin-bottom: 5px;
        }}
    </style>
</head>
<body>
    <div class="biorxiv-header">
        <div class="biorxiv-logo">bioRxiv</div>
        <div class="biorxiv-sub">PREPRINT | SYNTHETIC BIOLOGY & BIOENGINEERING</div>
    </div>
    {content}
</body>
</html>
"""

def prepare_html_content(md_text, lang="ru"):
    math_store = {}
    math_counter = 0

    # Preserve display math $$...$$
    def store_display_math(match):
        nonlocal math_counter
        key = f"MATHDISPLAYXYZ{math_counter}XYZ"
        math_counter += 1
        math_store[key] = match.group(0)
        return key

    # Preserve inline math $...$
    def store_inline_math(match):
        nonlocal math_counter
        key = f"MATHINLINEXYZ{math_counter}XYZ"
        math_counter += 1
        math_store[key] = match.group(0)
        return key

    md_text = re.sub(r"\$\$(.*?)\$\$", store_display_math, md_text, flags=re.DOTALL)
    md_text = re.sub(r"\$(.*?)\$", store_inline_math, md_text)

    # Format figure markdown ![Caption](path) into .figure-container
    def format_figure_md(match):
        alt_text = match.group(1).strip()
        rel_path = match.group(2).strip()
        abs_path = os.path.abspath(os.path.join(DOCS_DIR, rel_path))

        extra_class = ""
        if "graphical_abstract" in rel_path:
            extra_class = "graphical-abstract-img"

        caption_html = ""
        skip_captions = {"Graphical Abstract", "Графический абстракт"}
        if alt_text and alt_text not in skip_captions:
            caption_html = f'<div class="figure-caption"><strong>{alt_text}</strong></div>'

        return (
            f'<div class="figure-container">'
            f'<img class="{extra_class}" src="file://{abs_path}" alt="Figure">'
            f"{caption_html}</div>"
        )

    md_text = re.sub(r"!\[(.*?)\]\((.*?)\)", format_figure_md, md_text)

    # Convert Markdown to HTML
    html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "nl2br"])

    # Re-insert math blocks cleanly
    for key, math_code in math_store.items():
        html_body = html_body.replace(key, math_code)

    # Wrap metadata paragraph (RU / EN)
    html_body = re.sub(
        r"<p><strong>(Автор|Author):</strong>(.*?)</p>",
        r'<div class="meta-info"><strong>\1:</strong>\2</div>',
        html_body,
        flags=re.DOTALL,
    )
    # Collapse consecutive meta lines into one box when rendered as separate paragraphs
    html_body = re.sub(
        r'(<div class="meta-info">.*?</div>)(\s*<p><strong>(?:Контакты|Correspondence|Целев|Target|Дата|Date|Статус|Document|ПО|Software):</strong>.*?</p>)+',
        lambda m: m.group(0)
        .replace("</div>", "")
        .replace('<p><strong>', "<br><strong>")
        .replace("</p>", "")
        + "</div>",
        html_body,
        count=1,
        flags=re.DOTALL,
    )

    # Wrap Abstract / Аннотация
    html_body = re.sub(
        r"<h2>(Аннотация|Abstract)</h2>\s*<p>(.*?)</p>",
        r'<div class="abstract-box"><h2>\1</h2><p>\2</p></div>',
        html_body,
        flags=re.DOTALL,
    )

    # Basic table styling for readability in PDF
    html_body = html_body.replace(
        "<table>",
        '<table style="width:100%; border-collapse:collapse; font-size:9pt; margin:10px 0;">',
    )
    html_body = html_body.replace(
        "<th>",
        '<th style="border:1px solid #ccc; padding:4px 6px; background:#f6f8fa; text-align:left;">',
    )
    html_body = html_body.replace(
        "<td>",
        '<td style="border:1px solid #ddd; padding:4px 6px; vertical-align:top;">',
    )

    template = HTML_TEMPLATE
    if lang == "en":
        template = template.replace('lang="ru"', 'lang="en"')
        template = template.replace(
            "PREPRINT | SYNTHETIC BIOLOGY & BIOENGINEERING",
            "PREPRINT | SYNTHETIC BIOLOGY / BIOENGINEERING | v1.1",
        )

    return template.format(content=html_body)


def build_and_render_manuscript(lang="ru", screenshot_prefix=None):
    """Build HTML + PDF (+ page PNGs) for ru or en manuscript."""
    os.makedirs(REPORTS_DIR, exist_ok=True)
    lang = lang.lower()
    if lang not in ("ru", "en"):
        raise ValueError("lang must be 'ru' or 'en'")

    stem = f"manuscript_biorxiv_{lang}"
    md_file = os.path.join(DOCS_DIR, f"{stem}.md")
    html_file = os.path.join(DOCS_DIR, f"{stem}.html")
    pdf_file = os.path.join(DOCS_DIR, f"{stem}.pdf")
    prefix = screenshot_prefix or f"manuscript_{lang}_page_"

    if not os.path.isfile(md_file):
        raise FileNotFoundError(md_file)

    print(f"=== [{lang.upper()}] Step 1: Markdown → HTML ===")
    with open(md_file, "r", encoding="utf-8") as f:
        md_text = f.read()

    html_content = prepare_html_content(md_text, lang=lang)
    with open(html_file, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"  [OK] HTML: {html_file}")

    print(f"\n=== [{lang.upper()}] Step 2: Playwright → PDF ===")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1200, "height": 1600})
        page.goto(f"file://{html_file}")
        page.wait_for_timeout(3500)
        page.pdf(
            path=pdf_file,
            format="A4",
            print_background=True,
            margin={"top": "12mm", "bottom": "12mm", "left": "12mm", "right": "12mm"},
        )
        browser.close()
    size_kb = os.path.getsize(pdf_file) / 1024
    print(f"  [OK] PDF: {pdf_file} ({size_kb:.0f} KB)")

    print(f"\n=== [{lang.upper()}] Step 3: Page screenshots ===")
    doc = fitz.open(pdf_file)
    n_pages = doc.page_count
    screenshot_paths = []
    for page_num in range(n_pages):
        page = doc[page_num]
        pix = page.get_pixmap(dpi=150)
        img_name = f"{prefix}{page_num + 1}.png"
        img_path = os.path.join(REPORTS_DIR, img_name)
        pix.save(img_path)
        screenshot_paths.append(img_path)
        print(f"  [OK] {img_name}")
    doc.close()
    print(f"\n[SUCCESS] {lang.upper()}: {n_pages} pages → {pdf_file}")
    return {
        "pdf": pdf_file,
        "html": html_file,
        "pages": n_pages,
        "screenshots": screenshot_paths,
    }


if __name__ == "__main__":
    langs = sys.argv[1:] if len(sys.argv) > 1 else ["ru", "en"]
    results = []
    for lang in langs:
        results.append(build_and_render_manuscript(lang=lang))
    print("\n=== BUILD SUMMARY ===")
    for r in results:
        print(f"  {r['pdf']}  ({r['pages']} pages)")
