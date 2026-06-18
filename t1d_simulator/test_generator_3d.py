import sys
import os
from generator_3d import TPMSGenerator

def run_mesh_tests():
    print("=== Запуск автоматических тестов генератора 3D TPMS-структур ===")
    
    # Тест 1: Валидация сетки Gyroid на водонепроницаемость (watertight)
    print("Тест 1: Генерация гироида и проверка watertight-свойства...")
    gen_gyroid = TPMSGenerator(
        tpms_type="gyroid",
        box_size=600.0,
        unit_cell=200.0,
        thickness=0.0,
        resolution=60
    )
    mesh_g = gen_gyroid.build_mesh()
    
    assert mesh_g is not None, "Сетка Gyroid должна быть успешно построена."
    assert mesh_g.is_watertight, "Сетка Gyroid должна быть полностью замкнутой (watertight manifold)!"
    print(f"  [OK] Гироид сгенерирован: Вершин={len(mesh_g.vertices)}, Граней={len(mesh_g.faces)}")
    print(f"  [OK] Пористость гироида: {gen_gyroid.porosity:.2f}% | Удельная площадь SA/V: {gen_gyroid.sav_ratio:.2f} см⁻¹")
    
    # Тест 2: Монотонная зависимость пористости от параметра смещения (thickness)
    print("\nТест 2: Проверка монотонности зависимости пористости от смещения...")
    gen_thin = TPMSGenerator(tpms_type="gyroid", box_size=600.0, unit_cell=200.0, thickness=0.5, resolution=60)
    mesh_thin = gen_thin.build_mesh()
    
    gen_thick = TPMSGenerator(tpms_type="gyroid", box_size=600.0, unit_cell=200.0, thickness=-0.5, resolution=60)
    mesh_thick = gen_thick.build_mesh()
    
    print(f"  - Смещение t =  0.5 (тонкие стенки): Пористость = {gen_thin.porosity:.2f}%")
    print(f"  - Смещение t =  0.0 (стандарт):     Пористость = {gen_gyroid.porosity:.2f}%")
    print(f"  - Смещение t = -0.5 (толстые стенки): Пористость = {gen_thick.porosity:.2f}%")
    
    assert gen_thin.porosity > gen_gyroid.porosity, "Пористость при t=0.5 должна быть больше пористости при t=0.0."
    assert gen_gyroid.porosity > gen_thick.porosity, "Пористость при t=0.0 должна быть больше пористости при t=-0.5."
    print("  [OK] Закон монотонности изменения объема подтвержден.")
    
    # Тест 3: Генерация и экспорт тестовых STL
    print("\nТест 3: Экспорт тестовых файлов STL...")
    output_dir = "runs"
    os.makedirs(output_dir, exist_ok=True)
    
    gyroid_path = os.path.join(output_dir, "test_gyroid.stl")
    gen_gyroid.export_stl(gyroid_path)
    assert os.path.exists(gyroid_path), "Файл test_gyroid.stl должен существовать."
    
    schwarz_path = os.path.join(output_dir, "test_schwarz.stl")
    gen_schwarz = TPMSGenerator(tpms_type="schwarz_p", box_size=600.0, unit_cell=200.0, thickness=0.0, resolution=60)
    mesh_s = gen_schwarz.build_mesh()
    gen_schwarz.export_stl(schwarz_path)
    assert os.path.exists(schwarz_path), "Файл test_schwarz.stl должен существовать."
    print("  [OK] Файлы STL успешно экспортированы.")
    
    # Сравнение SA/V с плоским листом эквивалентного размера
    # Плоский лист размером 600x600x600 мкм при толщине стенки 600 мкм (сплошной куб):
    # SA/V куба = 6 / R = 6 / 0.03 cm = 200 см⁻¹
    # Но для листа толщиной 300 мкм (L=150) SA/V = 1 / L = 1 / 0.015 cm = 66.7 см⁻¹
    print(f"\n  Сравнение SA/V:")
    print(f"  - Гироид: {gen_gyroid.sav_ratio:.1f} см⁻¹")
    print(f"  - Поверхность Шварца P: {gen_schwarz.sav_ratio:.1f} см⁻¹")
    assert gen_gyroid.sav_ratio > 150.0, f"Удельная площадь гироида должна быть высокой, получено: {gen_gyroid.sav_ratio:.1f} см⁻¹"
    
    print("\n=== Все тесты генератора 3D-сеток успешно пройдены! ===")
    return True

if __name__ == "__main__":
    try:
        run_mesh_tests()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n  [ERROR] Тест провалена: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n  [ERROR] Непредвиденная ошибка: {e}", file=sys.stderr)
        sys.exit(1)
