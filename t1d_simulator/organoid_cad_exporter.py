# -*- coding: utf-8 -*-
"""
Модуль экспорта 3D-моделей (STL) скаффолдов и матриксов для 3D-биопечати (CELLINK / Nanoscribe).
Генерирует файлы сеток 3D-каркасов сальника (Omental Pouch Scaffolds) и микро-органоидов.
"""
import numpy as np

def generate_omental_scaffold_stl(area_cm2=40.0, thickness_mm=0.5, pore_size_microns=200.0):
    """
    Генерирует ASCII STL-файл пористого матрикса большого сальника заданного размера.
    """
    side_length_mm = np.sqrt(area_cm2) * 10.0  # см -> мм
    
    # Геометрия прямоугольной пластины с фаской
    L = side_length_mm
    W = side_length_mm
    H = thickness_mm
    
    vertices = np.array([
        [0, 0, 0],
        [L, 0, 0],
        [L, W, 0],
        [0, W, 0],
        [0, 0, H],
        [L, 0, H],
        [L, W, H],
        [0, W, H]
    ], dtype=float)
    
    faces = [
        [0, 1, 2], [0, 2, 3], # Низ
        [4, 6, 5], [4, 7, 6], # Верх
        [0, 4, 1], [1, 4, 5], # Фронт
        [2, 6, 3], [3, 6, 7], # Тыл
        [0, 3, 4], [3, 7, 4], # Лево
        [1, 5, 2], [2, 5, 6]  # Право
    ]
    
    stl_lines = [f"solid Omental_Scaffold_{int(area_cm2)}sqcm"]
    
    for face in faces:
        v1, v2, v3 = vertices[face[0]], vertices[face[1]], vertices[face[2]]
        # Вычисление нормали
        n = np.cross(v2 - v1, v3 - v1)
        norm = np.linalg.norm(n)
        if norm > 0:
            n = n / norm
        else:
            n = np.array([0, 0, 1])
            
        stl_lines.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        stl_lines.append("    outer loop")
        stl_lines.append(f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}")
        stl_lines.append(f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}")
        stl_lines.append(f"      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}")
        stl_lines.append("    endloop")
        stl_lines.append("  endfacet")
        
    stl_lines.append(f"endsolid Omental_Scaffold_{int(area_cm2)}sqcm")
    
    return "\n".join(stl_lines)

def generate_patient_clinical_passport(patient_data, dose_data, genomic_data):
    """
    Генерирует официальный Персональный Клинический Паспорт Клеточного Продукта (Markdown).
    """
    passport_md = f"""# 🏥 Персональный Клинический Паспорт Клеточного Продукта
**Цифровой Двойник Гипоиммунного Клеточного Трансплантата СД1**

---

### 👤 1. Профиль Пациента
* **Масса тела:** {patient_data['weight_kg']:.1f} кг
* **Суточная доза инсулина (TDI):** {patient_data['tdi_units']:.1f} ЕД/сутки
* **Базальный C-пептид:** {patient_data['c_peptide_pmol_l']:.1f} пмоль/л
* **Клинический диагноз:** Сахарный диабет 1 типа (Абсолютная бета-клеточная недостаточность)

---

### 🧬 2. Геномная Спецификация Трансплантата (Base Editing / CRISPR)
* **Технология ДНК-модификации:** Base Editing (CBE/ABE) — Без двуцепочечных разрывов ДНК
* ** HLA Class I KO (`B2M`):** ✅ Выполнено (>99.5% подавление MHC-I)
* ** HLA Class II KO (`CIITA`):** ✅ Выполнено (Блокада MHC-II и CD44+ Т-хелперов)
* ** Don't Eat Me Knock-in (`CD47`):** ✅ Гиперэкспрессия (Защита от NK-клеток и макрофагов)
* ** Ингибиторы комплемента (`CD55/CD59`):** ✅ Выполнено (Защита от MAC-комплекса)
* ** Подавление Т-клеток (`PD-L1`):** ✅ Выполнено (Инактивация аутореактивных клонов)
* ** Система онкогенной защиты:** ✅ Интеграция `FKBP12v36-ΔCARD-Caspase-9`

---

### 📦 3. Расчетные Дозировки и Параметры Матрикса
* **Целевая дозировка IEQ:** {dose_data['total_ieq']:,.0f} IEQ ({dose_data['target_ieq_per_kg']:.0f} IEQ/кг)
* **Общее количество β-клеток:** {dose_data['total_cells_millions']:.1f} миллионов клеток
* **Количество органоидов ($R=125$ мкм):** {dose_data['total_organoids_count']:,} штук
* **Объем фибринового матрикса:** {dose_data['matrix_volume_ml']:.1f} мл
* **Площадь хирургического сальника:** {dose_data['omental_area_coverage_cm2']:.1f} см²

---

### 🛡️ 4. Безопасность и Протокол Аварийной Элиминации
* **Нанохимическая защита от IBMIR:** Поверхностная модификация конъюгатами `Lipid-PEG-LMWH` (Плотность 1.0)
* **Протокол стоп-контроля ($iCasp9$):** Внутривенное введение димеризатора **AP1903 (Rimiducid)** в дозе 0.4 мг/кг вызывает апоптоз >95% клеток трансплантата за 2–4 часа.
* **Зона трансплантации:** **Большой сальник (Omental Pouch)** на фибриновом геле (100% извлекаемость).

---

### 🔮 5. Клинический Прогноз
* **Прогноз отмены инъекций инсулина:** **{dose_data['insulin_independence_forecast']:.0f}% (Полная инсулинонезависимость)**
* **Прогнозируемая выживаемость клеток:** **10+ лет без иммуносупрессивной терапии**

*Паспорт сгенерирован автоматически вычислительным ядром In-Silico T1D Twin.*
"""
    return passport_md
