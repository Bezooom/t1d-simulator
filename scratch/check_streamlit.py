import time
import os
from playwright.sync_api import sync_playwright

def main():
    print("Launching playwright...")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1080})
        print("Navigating...")
        page.goto("http://localhost:8501")
        print("Waiting 10 seconds for streamlit to load...")
        time.sleep(10)
        
        # Print all visible text
        print("Page text:")
        print(page.locator("body").inner_text())
        
        # Take a screenshot to inspect
        os.makedirs("scratch", exist_ok=True)
        page.screenshot(path="scratch/app_load_debug.png")
        print("Saved debug screenshot to scratch/app_load_debug.png")
        
        # Print H1 tags if any
        h1s = page.locator("h1").all()
        print(f"Found {len(h1s)} H1 tags:")
        for idx, h1 in enumerate(h1s):
            print(f"  {idx}: {h1.inner_text()}")
            
        browser.close()

if __name__ == "__main__":
    main()
