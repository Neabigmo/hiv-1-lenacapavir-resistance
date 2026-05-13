#!/usr/bin/env python3
"""
Convert HTML figures to high-quality PNG images using Playwright
"""

import asyncio
from pathlib import Path
from playwright.async_api import async_playwright

BASE_DIR = Path(__file__).parent.parent.parent
HTML_DIR = BASE_DIR / "scripts" / "revision_v2" / "html_figures"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures" / "html_renders"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

async def render_html_to_png(html_file, output_file, width=1200, height=1400):
    """Render HTML file to PNG using Playwright"""
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        page = await browser.new_page(viewport={'width': width, 'height': height})

        # Load HTML file
        file_url = html_file.absolute().as_posix().replace('\\', '/')
        await page.goto(f'file:///{file_url}')

        # Wait for rendering
        await page.wait_for_timeout(2000)

        # Take screenshot with high quality
        await page.screenshot(path=str(output_file), full_page=False, scale='device')

        await browser.close()
        print(f"[OK] Rendered: {output_file.name}")

async def main():
    """Render all HTML figures"""
    print("="*60)
    print("Rendering HTML Figures to PNG")
    print("="*60)

    figures = [
        ("figure4_structure_mechanisms.html", "figure4_ab.png", 1400, 700),
        ("figure5_resistance_pathway.html", "figure5_d.png", 900, 700),
        ("figure6_framework_bc.html", "figure6_bc.png", 1400, 700),
    ]

    for html_name, png_name, width, height in figures:
        html_file = HTML_DIR / html_name
        output_file = OUTPUT_DIR / png_name

        if html_file.exists():
            print(f"\nRendering {html_name}...")
            await render_html_to_png(html_file, output_file, width, height)
        else:
            print(f"[SKIP] {html_name} not found")

    print("\n" + "="*60)
    print("All HTML figures rendered!")
    print("="*60)

if __name__ == "__main__":
    asyncio.run(main())
