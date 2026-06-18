import torch
import numpy as np
import sys
import os

# Добавляем t1d_simulator в пути поиска модулей, если необходимо
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from gnn_pipeline import build_training_dataset, train_gnn_model, smiles_to_graph, BiocompatibilityGNN

def evaluate_molecules():
    print("=== Фаза 3: Обучение GNN и скрининг антифиброзных покрытий ===")
    
    # 1. Загрузка обучающего датасета
    dataset = build_training_dataset()
    print(f"Размер обучающей выборки: {len(dataset)} мономеров.")
    
    # 2. Обучение GNN
    epochs = 150
    print(f"Обучение GNN модели ({epochs} эпох)...")
    model, losses = train_gnn_model(dataset, epochs=epochs, lr=0.01, batch_size=8)
    
    # Верификация снижения ошибки
    initial_loss = losses[0]
    final_loss = losses[-1]
    loss_reduction = (initial_loss - final_loss) / initial_loss * 100
    print(f"Начальная ошибка (MSE): {initial_loss:.5f}")
    print(f"Конечная ошибка (MSE): {final_loss:.5f}")
    print(f"Снижение ошибки: {loss_reduction:.2f}%")
    
    assert loss_reduction > 70.0, f"Ошибка снизилась только на {loss_reduction:.2f}%, ожидалось >70%"
    print("  [OK] Модель успешно сошлась (снижение ошибки > 70%).")
    
    # 3. База данных кандидатов для скрининга (10 молекул)
    candidates = [
        {"name": "Zwitterion (Sulfobetaine SBAA)", "smiles": "C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]", "class": "Zwitterionic"},
        {"name": "Zwitterion (Carboxybetaine CBAA)", "smiles": "C=CC(=O)NCC[N+](C)(C)CC(=O)[O-]", "class": "Zwitterionic"},
        {"name": "PEG8-acrylate", "smiles": "C=CC(=O)OCCOCCOCCOCCOCCOCCOCCOCCO C", "class": "PEGylated"},
        {"name": "PEG4-methacrylate", "smiles": "CC(=C)C(=O)OCCOCCOCCOCCO C", "class": "PEGylated"},
        {"name": "N,N-dimethylacrylamide (DMAA)", "smiles": "C=CC(=O)N(C)C", "class": "Neutral Hydrophilic"},
        {"name": "Hydroxyethyl methacrylate (HEMA)", "smiles": "CC(=C)C(=O)OCCO", "class": "Neutral Hydrophilic"},
        {"name": "Methyl methacrylate (MMA)", "smiles": "CC(=C)C(=O)OC", "class": "Hydrophobic"},
        {"name": "Styrene", "smiles": "C=CC1=CC=CC=C1", "class": "Hydrophobic"},
        {"name": "METAC (quaternary cationic)", "smiles": "CC(=C)C(=O)OCC[N+](C)(C)C", "class": "Cationic"},
        {"name": "Pentafluoropropyl methacrylate (PFPMA)", "smiles": "CC(=C)C(=O)OCC(F)(F)C(F)(F)F", "class": "Fluorinated Hydrophobic"}
    ]
    
    # Оценка кандидатов
    model.eval()
    results = []
    max_fibrosis_thickness = 150.0 # мкм
    
    with torch.no_grad():
        for cand in candidates:
            data = smiles_to_graph(cand["smiles"])
            if data is None:
                print(f"  [ERROR] Не удалось распарсить SMILES для {cand['name']}")
                continue
            
            # Создаем тензор batch для одного графа
            batch = torch.zeros(data.x.size(0), dtype=torch.long)
            pred = model(data.x, data.edge_index, batch)
            biocompatibility = float(pred[0, 0].item())
            
            # Эмпирическое отображение на толщину фиброза L_fib (мкм)
            L_fib = max(0.0, max_fibrosis_thickness * (1.0 - biocompatibility))
            
            results.append({
                "name": cand["name"],
                "smiles": cand["smiles"],
                "class": cand["class"],
                "biocompatibility": biocompatibility,
                "L_fib": L_fib
            })
            
    # Ранжирование по убыванию биосовместимости
    results.sort(key=lambda x: x["biocompatibility"], reverse=True)
    
    print("\n=== Результаты скрининга кандидатов ===")
    print(f"{'Название покрытия':<35} | {'Класс':<20} | {'Индекс совместимости':<21} | {'Ожидаемый фиброз L_fib (мкм)':<28} |")
    print("-" * 115)
    for r in results:
        print(f"{r['name']:<35} | {r['class']:<20} | {r['biocompatibility']:<21.4f} | {r['L_fib']:<28.1f} |")
        
    print("\n=== Топ-3 антифиброзных покрытия ===")
    for i, r in enumerate(results[:3]):
        print(f"{i+1}. {r['name']} ({r['class']}) - Совместимость: {r['biocompatibility']:.4f}, Фиброз: {r['L_fib']:.1f} мкм")
        
    # Проверка биологического ранжирования
    zwitter_scores = [r["biocompatibility"] for r in results if r["class"] == "Zwitterionic"]
    peg_scores = [r["biocompatibility"] for r in results if r["class"] == "PEGylated"]
    hydrophobic_scores = [r["biocompatibility"] for r in results if r["class"] == "Hydrophobic"]
    cationic_scores = [r["biocompatibility"] for r in results if r["class"] == "Cationic"]
    
    assert min(zwitter_scores) > max(peg_scores), "Ошибка ранжирования: Цвиттер-ионы должны быть совместимее ПЭГ"
    assert min(peg_scores) > max(hydrophobic_scores), "Ошибка ранжирования: ПЭГ должен быть совместимее гидрофобных мономеров"
    assert min(hydrophobic_scores) > min(cationic_scores), "Ошибка ранжирования: Гидрофобные должны быть совместимее катионных"
    print("  [OK] Биологическое ранжирование подтверждено.")
    
    # Сохраняем веса обученной модели
    save_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "biocompatibility_gnn.pt")
    torch.save(model.state_dict(), save_path)
    print(f"  [OK] Веса обученной GNN сохранены в {save_path}")

if __name__ == "__main__":
    try:
        evaluate_molecules()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n  [ERROR] Тест провален: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n  [ERROR] Непредвиденная ошибка: {e}", file=sys.stderr)
        sys.exit(1)
