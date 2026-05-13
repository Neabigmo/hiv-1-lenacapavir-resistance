#!/usr/bin/env python3
"""
HTML to Image Converter
将HTML文件转换为高质量的PNG和PDF图片
"""

import os
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright
import time

def html_to_image(html_path, output_dir, format='both'):
    """
    将HTML文件转换为图片

    Args:
        html_path: HTML文件路径
        output_dir: 输出目录
        format: 'png', 'pdf', 或 'both'
    """
    html_path = Path(html_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_name = html_path.stem

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()

        # 加载HTML文件
        page.goto(f'file:///{html_path.absolute().as_posix()}')

        # 等待页面完全加载
        page.wait_for_load_state('networkidle')
        time.sleep(1)  # 额外等待确保渲染完成

        # 生成PNG
        if format in ['png', 'both']:
            png_path = output_dir / f'{base_name}.png'
            page.screenshot(
                path=str(png_path),
                full_page=True,
                type='png'
            )
            print(f'[OK] Generated PNG: {png_path}')

        # 生成PDF
        if format in ['pdf', 'both']:
            pdf_path = output_dir / f'{base_name}.pdf'
            page.pdf(
                path=str(pdf_path),
                format='A4',
                print_background=True,
                margin={'top': '0', 'right': '0', 'bottom': '0', 'left': '0'}
            )
            print(f'[OK] Generated PDF: {pdf_path}')

        browser.close()

def main():
    # 设置路径
    script_dir = Path(__file__).parent
    html_dir = script_dir / 'html_figures'
    output_dir = Path('H:/2026try/4.20/manuscript/figures/revision_v2/html_renders')

    # 获取所有HTML文件
    html_files = sorted(html_dir.glob('*.html'))

    if not html_files:
        print(f'No HTML files found in {html_dir}')
        return

    print(f'Found {len(html_files)} HTML files to convert')
    print(f'Output directory: {output_dir}\n')

    # 转换每个HTML文件
    for html_file in html_files:
        print(f'Converting: {html_file.name}')
        try:
            html_to_image(html_file, output_dir, format='both')
        except Exception as e:
            print(f'[ERROR] Error converting {html_file.name}: {e}')

    print(f'\n[OK] Conversion complete! Check {output_dir}')

if __name__ == '__main__':
    main()
