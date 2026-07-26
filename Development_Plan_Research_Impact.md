# План развития: исследования и вклад в терапию СД1

**Проект:** In Silico β-cell / hypoimmune organoid digital twin  
**Статус базы:** 10 фаз аудита (PDE/PINN, TPMS, GNN, OGM, ангиогенез, Phase 10 Biomimesis), Streamlit-симулятор, научный бриф  
**Цель документа:** зафиксировать, *как копать глубже* и *как реально помогать* лечению сахарного диабета 1 типа, без иллюзии «готового клеточного продукта»  
**Дата:** 2026-07-17  
**Актуальный план исполнения:** [`ROADMAP.md`](ROADMAP.md) **v3.0** (2026-07-26) — M0–M3 in-repo закрыты; внешние шаги + M4–M5. Этот файл — стратегический обзор impact-треков.

---

## 1. Исходная позиция

### 1.1. Что уже есть

| Артефакт | Содержание |
|----------|------------|
| Цифровой двойник | `t1d_simulator/`: O₂ PDE/PINN, TPMS Gyroid/Schwarz, GNN-антифиброз, VEGF, organoid (Phase 10) |
| Научные документы | Scientific Brief, Advanced Bio Design Solutions, Critical Analysis (фазы 6–10), Development Roadmap |
| Ключевой вывод | Макроинкапсуляция системно ограничена (Krogh, FBR, MWCO/IgG, death window, OGM/ROS); финальная парадигма — hypoimmune organoids + васкуляризация + suicide switch + anti-IBMIR + omentum |

### 1.2. Чего нет (и это нормально на текущем этапе)

- Wet-lab данных (клетки, животные, клиника)
- Независимого peer-review / preprint
- IP на линию iPSC или edit-set
- Калибровки параметров модели по экспериментальным paper-benchmarks
- Публичного demo / EN-упаковки для outreach

### 1.3. Принцип рычага

| Уровень | Влияние на лечение | Нужна ли своя lab |
|---------|-------------------|-------------------|
| **A.** Лучшие модели + open tools | Ускоряет чужие дизайны, режет тупики | Нет |
| **B.** Валидация на публичных данных | Доверие, цитируемость, партнёры | Нет |
| **C.** Co-authorship с wet-lab | Прямой вклад в протокол | Партнёр — да, мы — нет |
| **D.** Свой клеточный продукт | Максимум, годы и капитал | Да |

**Оптимальный путь малой группы: A → B → C.**  
Фаза 5 из `Development_Roadmap.md` (heparin coat, iCasp9 in vitro, WGS) — лабораторный трек; без партнёра туда не идём в одиночку.

---

## 2. Стратегические треки исследований

### Трек A — Сделать twin правдивым (высший ROI на текущей базе)

Не «ещё одна концептуальная фаза», а **калибровка, недостающие failure modes, decision tables**.

| # | Задача | Зачем | Выход |
|---|--------|-------|-------|
| A1 | **Literature-calibrated parameters** | Сейчас сила — уравнения; слабость — привязка к эксперименту | Таблица \(V_{max}\), \(K_M\), \(D_{O_2}\), \(L_{fib}\), time-to-vascularization (mouse/human), IBMIR loss % из ≥20 papers |
| A2 | **Benchmark reproduction** | «Воспроизводим ли Fig. 3 из paper X?» | Отчёт: predicted vs reported viability / pO₂ |
| A3 | **IBMIR 0–48 h module** | Главный клинический killer «голых» островков; в коде organoid-часть пока тонкая | ODE: TF → thrombin → clot → O₂ drop → 24–48 h viability |
| A4 | **Site comparison** | Portal vs omentum vs SQ — язык трансляции | Decision table: survival, steatosis risk, surgical complexity |
| A5 | **Virtual dose / mass** | Связь islet mass ↔ insulin independence (упрощённо) | Сценарии «therapeutic dose» с честными assumptions |
| A6 | **Public omics** | scRNA-seq SC-islets / native islets (GEO) | Hypoxia + immune genes; input для edit-set narrative |
| A7 | **In silico edit / off-target** | Язык CRISPR Tx / Sana | Minimal edit set (B2M, CIITA, CD47, CD55, iCasp9) + risk scores |
| A8 | **Code quality** | Без этого open science и партнёры невозможны | `parameters.yaml`, tests, EN README, reproducible scripts |

**Приоритет Q1:** A1 → A2 → A3 → A4 → A8.

### Трек B — Не только organoid (near-term польза пациентам)

Клеточный cure — горизонт 5–15 лет. Параллельный портфель:

| Направление | Почему важно | Skill fit |
|-------------|--------------|-----------|
| AID / closed-loop algorithms | Уже в клинике | ML, control, sim |
| CGM prediction / hypo prevention | Снижает риск и страх | time-series ML |
| Drug repurposing / immunotherapy meta-analysis | Newly diagnosed | NLP + stats |
| Patient stratification (C-peptide, honeymoon) | Кого лечить чем | data science |

**Правило:** long-shot = organoid twin; near-term = AID/CGM/data. Оба законны.

### Трек C — Открытая наука (мультипликатор)

| # | Задача | Выход |
|---|--------|-------|
| C1 | Preprint (EN) | bioRxiv: *Multiphysics failure modes of β-cell encapsulation: from Krogh limit to hypoimmune organoids* |
| C2 | Public repo | GitHub: clean simulator + DOI (Zenodo) |
| C3 | Interactive demo | Streamlit Cloud / Hugging Face: death window, OGM, organoid modes |
| C4 | Methods + limitations | Честный critical analysis → methods section (gold for trust) |

Без C1–C3 cold outreach почти мёртв.

### Трек D — Выход к wet-lab и индустрии

| # | Задача | Как |
|---|--------|-----|
| D1 | Academic partners | PubMed: islet encapsulation / hypoimmune iPSC / omental pouch; offer: free run на *их* geometry |
| D2 | Mid biotech (не только Vertex) | Sana, CRISPR Tx, Seraxis, Sernova, Sigilon-track — персональный 1-pager |
| D3 | Grants / challenges | Breakthrough T1D, uni interdisciplinary IT+biophysics |
| D4 | Affiliation / real collab | Не только шапка в PDF — реальные письма в lab biophysics / bioengineering |

**Позиционирование:** не «у нас cure», а  
*independent multiphysics failure analysis + design-space screening; offering model co-development*.

---

## 3. 12-месячный план

### Q1 — Калибровка и инженерия модели

- [ ] Таблица параметров из ≥10–20 ключевых papers (A1)
- [ ] Минимум 2–3 benchmark reproduction (A2)
- [ ] Модуль IBMIR 0–48 h (A3)
- [ ] Portal vs omentum vs SQ comparison (A4)
- [ ] `parameters.yaml` + unit tests + EN README (A8)
- [ ] Черновик EN abstract preprint (C1 start)

**Артефакты:** `parameters/` или `data/literature_params.yaml`; `organoid_simulator` + IBMIR; internal validation report.

### Q2 — Open science

- [ ] Preprint на bioRxiv (C1)
- [ ] Публичный demo (C3)
- [ ] GitHub + Zenodo DOI (C2)
- [ ] Опционально: scRNA short note / edit-set side project (A6–A7)

**Артефакты:** preprint PDF, public URL demo, citeable code.

### Q3 — Outreach и near-term track

- [ ] 15 targeted писем: labs + 3–5 biotech (D1–D2)
- [ ] Цель: 2–3 содержательных ответа / scientific conversation
- [ ] Опционально: mini-project AID/CGM (Трек B)

**Артефакты:** EN 1-pager + 5-page brief; email log; (opt) CGM/AID notebook.

### Q4 — Углубление партнёрства

- [ ] Joint figure / co-authorship **или** second preprint
- [ ] Grant application (computational digital twin for cell therapy)
- [ ] Обновление модели по feedback партнёров

**Артефакты:** shared paper / grant / v2 model release.

---

## 4. Что *не* делать (анти-план)

1. Ещё 5 манифестов «Ultimate Biomimesis» без данных — diminishing returns.  
2. Притворяться wet-lab / clinical readiness.  
3. Только cold email в Vertex BD без preprint/demo.  
4. Ждать «полную 3D CFD» перед публикацией — лучше calibrated 1D/2D + honest limits.  
5. Мультиплекс CRISPR «протокол» без доступа к iPSC — останется текстом.  
6. Подменять помощь пациентам только мечтой о cure: near-term треки (AID/CGM) тоже считаются.

---

## 5. Outreach: кому и с чем

### 5.1. Приоритет адресатов

| Приоритет | Кому | Почему |
|-----------|------|--------|
| 1 | **Sana Biotechnology** (HIP / SC451) | Максимальный overlap с Phase 10 |
| 2 | **CRISPR Therapeutics** (CTX211) | Gene-edit / hypoimmune safety modelling |
| 3 | **Seraxis, Sernova, Sigilon/Lilly-track** | Ещё в encapsulation — ценны фазы 1–9 (failure modes) |
| 4 | **Vertex** | Лидер; только через intro / conference / после preprint |
| 5 | **Academic PI** (islet bioengineering, biomaterials) | Быстрее ответ; путь к co-authorship |

### 5.2. Пакет документов для отправки

1. **EN 1-page executive summary**  
2. **Scientific Brief (EN, 5–6 стр.)** — тон recommendations, не «предписания»  
3. **Live demo / screencast 5 min**  
4. (По запросу) NDA + code / full methods  

### 5.3. Формула письма (шаблон смысла)

> We built a multiphysics digital twin of β-cell replacement designs (O₂ transport, FBR, angiogenesis lag, OGM toxicity, IBMIR).  
> We quantitatively map failure modes of macroencapsulation and design countermeasures aligned with hypoimmune organoid approaches.  
> Offer: free screening of *your* geometry / density / site parameters; co-development of calibrated modules.

---

## 6. Метрики успеха (не vanity)

| Метрика | Целевой ориентир (12 мес.) |
|---------|---------------------------|
| Papers reproduced / calibrated | ≥3 |
| Preprint | 1 |
| Public demo uptime | да |
| Lab/biotech meaningful replies | ≥2 |
| Code tests on core solvers | green CI или хотя бы local test suite |
| Optional near-term patient-facing tool | 0 или 1 (не блокер) |

**Не метрики:** число MD-файлов, число «фаз», громкость формулировок.

---

## 7. Связь с существующими документами

| Документ | Роль |
|----------|------|
| **`ROADMAP.md`** | **Единая живая дорожная карта (вехи, кварталы, kickstart)** |
| `Development_Roadmap.md` | Исторический roadmap фаз 1–5 симулятора + prep to in vivo |
| `Scientific_Brief_Biomaterials_Lab.md` | Научное обоснование Phase 10 для lab/R&D |
| `Advanced_Bio_Design_Solutions.md` | Полный научный доклад (PDE, TPMS, GNN, OGM, …) |
| `Critical_Analysis_*.md` | Честные пределы модели — основа methods/limitations |
| **`Development_Plan_Research_Impact.md`** (этот файл) | **Стратегия развития: research + open science + impact + outreach** |

Обновлять этот план при закрытии квартальных чекбоксов (дата + short changelog внизу).

---

## 8. Недельный kickstart (минимальный старт)

| День / слот | Действие |
|-------------|----------|
| 1–2 | Завести `data/literature_params.yaml` + 5 papers по O₂ / FBR |
| 3–4 | Набросать IBMIR ODE stub в `organoid_simulator.py` |
| 5 | EN abstract (250–300 слов) |
| 6 | Список 5 lab emails + персональный hook |
| 7 | README EN: что модель может / не может |

---

## 9. Личный горизонт (если путь на годы)

1. **Скиллы:** computational biology + mass transport + immunology basics (база уже есть).  
2. **Якорь:** uni lab collab / PhD / internship в cell therapy или diabetes tech.  
3. **Ниша:** «переводчик» между biomaterials engineers и immunologists.  
4. **Два горизонта:**  
   - *near-term* — CGM/AID, education, open tools;  
   - *long-term* — digital twins for cell replacement.  

**Помогать ≠ обязательно вырастить β-клетку своими руками.**  
**Помогать = ускорить тех, кто выращивает, и облегчить жизнь тем, кто ждёт.**

---

## 10. Ближайшие опции реализации (когда готов код/текст)

| ID | Работа | Тип |
|----|--------|-----|
| **A** | Каркас `parameters.yaml` + literature calibration script | Code |
| **B** | Модуль IBMIR 0–48 h | Code |
| **C** | EN abstract + structure bioRxiv | Writing |
| **D** | Шаблоны писем 5 academic labs + 1-pager | Outreach |

Рекомендуемый порядок: **A → B → C → D**.

---

## Changelog

| Дата | Изменение |
|------|-----------|
| 2026-07-17 | Первая версия плана развития (research + impact + 12-month roadmap) |
| 2026-07-23 | Ссылка на единый `ROADMAP.md` как источник правды по вехам |
