import time
import os
import sys
from playwright.sync_api import sync_playwright

def run_automation():
    artifact_dir = "/home/bezoom/.gemini/antigravity/brain/44418490-fbde-44a6-a563-7f18b664eb0a"
    screenshot_1d_path = os.path.join(artifact_dir, "screenshot.png")
    screenshot_gnn_path = os.path.join(artifact_dir, "screenshot_gnn.png")
    screenshot_vegf_path = os.path.join(artifact_dir, "screenshot_vegf.png")
    
    print("Инициализация Playwright...")
    with sync_playwright() as p:
        print("Запуск браузера Chromium (headless)...")
        browser = p.chromium.launch(headless=True)
        
        page = browser.new_page(viewport={"width": 1440, "height": 1080})
        
        url = "http://localhost:8501"
        print(f"Переход по адресу: {url} ...")
        page.goto(url)
        
        print("Ожидание появления заголовка H1...")
        try:
            page.wait_for_selector('h1', timeout=20000)
            print("Приложение загружено! Текущий заголовок:", page.title())
        except Exception as e:
            print("Таймаут ожидания H1:", e)
            
        time.sleep(3) # Ждем стабилизации Plotly
        
        # Шаг 1: Тестируем GNN режим
        print("Переключение в режим ML-подбора покрытий (GNN)...")
        try:
            page.click('text="🧪 ML-подбор антифиброзных покрытий (GNN)"')
            print("Успешно переключено на вкладку GNN.")
            time.sleep(5) # Ждем загрузки GNN и отрисовки RDKit
            
            print(f"Сохранение снимка экрана GNN в {screenshot_gnn_path} ...")
            page.screenshot(path=screenshot_gnn_path, full_page=True)
            
            print("Экспорт L_fib в симулятор...")
            page.click('button:has-text("Применить это покрытие")')
            time.sleep(2)
            
        except Exception as e:
            print("Ошибка при автоматизации GNN вкладки:", e)
            
        # Шаг 2: Тестируем режим Неоваскуляризации (VEGF)
        print("Переключение в режим Неоваскуляризации (VEGF)...")
        try:
            page.click('text="🩸 Неоваскуляризация (VEGF / Ангиогенез)"')
            print("Успешно переключено на вкладку VEGF.")
            time.sleep(5) # Ждем расчета и отрисовки Plotly графиков
            
            print(f"Сохранение снимка экрана VEGF в {screenshot_vegf_path} ...")
            page.screenshot(path=screenshot_vegf_path, full_page=True)
            
        except Exception as e:
            print("Ошибка при автоматизации VEGF вкладки:", e)
            
        # Шаг 3: Возвращаемся в 1D режим и проверяем импортированное L_fib
        print("Возврат в режим 1D симуляции...")
        try:
            page.click('text="1D Симуляция диффузии O₂"')
            print("Успешно переключено на вкладку 1D.")
            time.sleep(3)
            
            print(f"Сохранение снимка экрана 1D симуляции в {screenshot_1d_path} ...")
            page.screenshot(path=screenshot_1d_path, full_page=True)
            
        except Exception as e:
            print("Ошибка при автоматизации 1D вкладки после экспорта:", e)
            
        # Шаг 4: Тестируем режим мини-органоидов (Фаза 10)
        print("Переключение в режим мини-органоидов (Фаза 10)...")
        try:
            page.click('text="🧫 Мини-органоиды (Фаза 10: Biomimesis)"')
            print("Успешно переключено на вкладку мини-органоидов.")
            time.sleep(5)
            
            screenshot_organoid_path = os.path.join(artifact_dir, "screenshot_organoid.png")
            print(f"Сохранение снимка экрана мини-органоидов в {screenshot_organoid_path} ...")
            page.screenshot(path=screenshot_organoid_path, full_page=True)
            
        except Exception as e:
            print("Ошибка при автоматизации вкладки мини-органоидов:", e)
            
        browser.close()
        print("Автоматизация успешно завершена.")

if __name__ == "__main__":
    run_automation()
