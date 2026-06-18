import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from torch_geometric.data import Data, DataLoader

# --- 1. Конвертер SMILES -> Граф Молекулы (RDKit + PyTorch Geometric) ---
def smiles_to_graph(smiles, y_val=None):
    """
    Парсит SMILES строку и преобразует ее в граф PyTorch Geometric Data.
    Node features: One-hot тип атома (C, N, O, S, P, F, Cl, Br, I) + заряд + ароматичность.
    Edge features: Bond type.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
        
    atom_list = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']
    node_feats = []
    
    # 1. Формирование признаков вершин (атомов)
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        # One-hot кодирование типа атома
        one_hot = [1.0 if symbol == s else 0.0 for s in atom_list]
        
        # Добавляем заряд и ароматичность
        charge = float(atom.GetFormalCharge())
        is_aromatic = 1.0 if atom.GetIsAromatic() else 0.0
        
        # Итоговый вектор признаков атома (размерность 11)
        feat = one_hot + [charge, is_aromatic]
        node_feats.append(feat)
        
    x = torch.tensor(node_feats, dtype=torch.float32)
    
    # 2. Формирование ребер (химических связей)
    edge_indices = []
    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        # Двунаправленные ребра (для неориентированного графа)
        edge_indices.append([start, end])
        edge_indices.append([end, start])
        
    edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    if len(edge_indices) == 0:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        
    # Целевой признак (индекс биосовместимости в [0, 1])
    y = torch.tensor([[y_val]], dtype=torch.float32) if y_val is not None else None
    
    return Data(x=x, edge_index=edge_index, y=y)

# --- 2. Архитектура графовой нейросети (GNN) ---
class BiocompatibilityGNN(nn.Module):
    """
    Графовая нейросеть для предсказания биосовместимости покрытия по структуре молекулы.
    """
    def __init__(self, input_dim=11, hidden_dim=32):
        super().__init__()
        from torch_geometric.nn import GCNConv, global_mean_pool
        
        # 3 слоя свертки GCN для message-passing
        self.conv1 = GCNConv(input_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, hidden_dim)
        
        # Полносвязный классификатор (MLP)
        self.fc1 = nn.Linear(hidden_dim, hidden_dim // 2)
        self.fc2 = nn.Linear(hidden_dim // 2, 1)
        
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(0.1)
        
    def forward(self, x, edge_index, batch):
        from torch_geometric.nn import global_mean_pool
        # 1. Message passing
        h = self.relu(self.conv1(x, edge_index))
        h = self.relu(self.conv2(h, edge_index))
        h = self.relu(self.conv3(h, edge_index))
        
        # 2. Глобальный пулинг графа (агрегируем вершины в один вектор)
        h_graph = global_mean_pool(h, batch)
        
        # 3. Классификация
        h_out = self.relu(self.fc1(h_graph))
        h_out = self.dropout(h_out)
        out = self.sigmoid(self.fc2(h_out))
        return out

# --- 3. Куратор обучающего датасета молекул ---
def build_training_dataset():
    """
    Генерирует расширенную базу из 52 референсных молекул/мономеров пяти классов
    с установленными индексами биосовместимости (сопротивления FBR).
    """
    raw_data = []
    
    # 1. Цвиттер-ионы (Индекс: 0.90 - 0.98) - максимальная биосовместимость
    zwitterions = [
        ("CC(=C)C(=O)OCCN(C)(C)CCCS(=O)(=O)[O-]", 0.96),  # SBMA
        ("C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]", 0.98),  # SBAA
        ("CC(=C)C(=O)OCCN(C)(C)CC(=O)[O-]", 0.95),       # CBMA
        ("C=CC(=O)NCC[N+](C)(C)CC(=O)[O-]", 0.97),       # CBAA
        ("CC(=C)C(=O)OCOP(=O)([O-])OCC[N+](C)(C)C", 0.94), # MPC
        # Аугментация структуры (вариации спейсеров и алкильных цепей)
        ("CC(=C)C(=O)OCCN(C)(C)CCCCS(=O)(=O)[O-]", 0.95),
        ("C=CC(=O)NCC[N+](C)(C)CCCCS(=O)(=O)[O-]", 0.97),
        ("CC(=C)C(=O)OCCN(CC)(CC)CCCS(=O)(=O)[O-]", 0.94),
        ("C=CC(=O)N1CCN(CCCS(=O)(=O)[O-])CC1", 0.96),
        ("C=CC(=O)OCCN(C)(C)CCCS(=O)(=O)[O-]", 0.96),
        ("CC(=C)C(=O)OCOP(=O)(OCC[N+](CC)(CC)C)[O-]", 0.92)
    ]
    raw_data.extend(zwitterions)
    
    # 2. ПЭГилированные мономеры (Индекс: 0.75 - 0.85) - высокая биосовместимость
    pegs = [
        ("CC(=C)C(=O)OCCOCCOCCOCCOCCO C", 0.83),  # OEGMA-5
        ("C=CC(=O)OCCOCCOCCOCCO C", 0.81),       # PEG4-acrylate
        ("C=CC(=O)OCCOCCOCCO C", 0.80),          # PEG3-acrylate
        ("CC(=C)C(=O)OCCOCCO C", 0.78),
        # Аугментированные цепочки PEG разной длины
        ("C=CC(=O)OCCOCCOCCOCCOCCO C", 0.82),
        ("C=CC(=O)OCCOCCOCCOCCOCCOCCO C", 0.83),
        ("C=CC(=O)OCCOCCOCCOCCOCCOCCOCCO C", 0.84),
        ("CC(=C)C(=O)OCCOCCOCCO C", 0.79),
        ("CC(=C)C(=O)OCCOCCOCCOCCO C", 0.81),
        ("CC(=C)C(=O)OCCOCCOCCOCCOCCOCCO C", 0.83)
    ]
    raw_data.extend(pegs)
    
    # 3. Нейтральные гидрофильные (Индекс: 0.50 - 0.65) - умеренная совместимость
    hydrophilic = [
        ("CC(=C)C(=O)OCCO", 0.62),      # HEMA
        ("C=CC(=O)N", 0.55),             # Acrylamide
        ("C=CC(=O)NC(C)C", 0.58),        # NIPAM
        ("C=CCN1CCCC1=O", 0.60),         # N-vinylpyrrolidone
        # Аугментация
        ("CC(=C)C(=O)NCC(O)CO", 0.64),
        ("C=CC(=O)N(C)C", 0.57),
        ("CC(=C)C(=O)NC(C)C", 0.59),
        ("C=CC(=O)N1CCOCC1", 0.61),
        ("C=CCN1CCOCC1=O", 0.58),
        ("CC(=C)C(=O)OCCN1CCOCC1", 0.63),
        ("C=CC(=O)OCCN(C)C", 0.54)
    ]
    raw_data.extend(hydrophilic)
    
    # 4. Гидрофобные мономеры (Индекс: 0.15 - 0.35) - низкая совместимость (FBR)
    hydrophobic = [
        ("CC(=C)C(=O)OC", 0.28),         # MMA
        ("C=CC1=CC=CC=C1", 0.22),        # Styrene
        ("CC(=C)C(=O)OCCCC", 0.24),      # BMA
        ("CC(=C)C(=O)OCC(F)(F)F", 0.18),  # TFEMA
        ("CC(=C)C(=O)OC(C(F)(F)F)C(F)(F)F", 0.15),
        # Аугментация
        ("C=CC(=O)OC", 0.26),
        ("C=CC(=O)OCCCC", 0.23),
        ("C=CC(=O)OCCCCCCC", 0.20),
        ("CC(=C)C(=O)OCCCCCCCCCCCC", 0.16), # LMA
        ("CC(=C)C(=O)OCCCCCC", 0.22),
        ("C=CC1=CC=C(C)C=C1", 0.21)
    ]
    raw_data.extend(hydrophobic)
    
    # 5. Катионные мономеры (Индекс: 0.05 - 0.15) - крайне низкая (адгезия, цитотоксичность)
    cationic = [
        ("CC(=C)C(=O)OCC[N+](C)(C)C", 0.12),  # METAC
        ("C=CC(=O)OCC[N+](C)(C)C", 0.10),     # AETAC
        ("C=C[N+](C)(C)CC=C", 0.08),          # DADMAC
        # Аугментация
        ("CC(=C)C(=O)OCCN(C)(C)C", 0.13),
        ("CC(=C)C(=O)OCC[N+](CC)(CC)CC", 0.09),
        ("C=CC(=O)OCC[N+](CC)(CC)CC", 0.08),
        ("C=CC(=O)NCC[N+](C)(C)C", 0.11),
        ("CC(=C)C(=O)NCC[N+](C)(C)C", 0.12),
        ("C=CCN(C)(C)C", 0.07),
        ("CC(=C)C(=O)OCC[N+](C)(C)CCCCCC", 0.06)
    ]
    raw_data.extend(cationic)
    
    # Преобразуем SMILES в графы PyG
    dataset = []
    for smiles, y_val in raw_data:
        data = smiles_to_graph(smiles, y_val)
        if data is not None:
            dataset.append(data)
            
    return dataset

# --- 4. Функция обучения GNN ---
def train_gnn_model(dataset, epochs=150, lr=0.01, batch_size=8):
    """
    Обучает графовую нейросеть предсказывать биосовместимость полимера.
    """
    torch.manual_seed(42)
    
    # Инициализация загрузчика данных
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    model = BiocompatibilityGNN()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    criterion = nn.MSELoss()
    
    # Записываем историю ошибок
    losses = []
    
    model.train()
    for epoch in range(epochs):
        epoch_loss = 0.0
        for batch in loader:
            optimizer.zero_grad()
            out = model(batch.x, batch.edge_index, batch.batch)
            loss = criterion(out, batch.y)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * batch.num_graphs
            
        epoch_loss /= len(dataset)
        losses.append(epoch_loss)
        
    return model, losses
