import os
import yaml

DEFAULT_PARAMS_PATH = os.path.join(os.path.dirname(__file__), "parameters.yaml")
LITERATURE_PARAMS_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "literature_params.yaml")

DEFAULT_FALLBACK_PARAMS = {
    "solubility": 1.34e-9,
    "v_max": 1.2e-16,
    "v_max_human": 1.5e-16,
    "v_max_macrophage": 3.0e-16,
    "k_m": 0.5,
    "k_m_insulin": 5.0,
    "hydrogels": {
        "water": {"name": "Чистая вода (контроль)", "D": 3.0e-5},
        "alginate_1%": {"name": "1% Альгинат натрия", "D": 2.1e-5},
        "alginate_2%": {"name": "2% Альгинат натрия (стандарт)", "D": 1.5e-5},
        "fibrosis_tissue": {"name": "Плотный фиброз", "D": 1.0e-5}
    },
    "implantation_sites": {
        "arterial": {"name": "Артериальное русло (прямая перфузия / AV-петля)", "pO2": 95.0, "description": "Идеальные условия, высокая оксигенация."},
        "omental_pouch": {"name": "Сальник (Omental pouch, реваскуляризованный)", "pO2": 55.0, "description": "Богатая сосудистая сеть, доступность извлечения, ниже риск IBMIR."},
        "venous": {"name": "Венозная сеть (воротная вена / печень)", "pO2": 40.0, "description": "Умеренная оксигенация, прямой контакт с кровью (высокий IBMIR)."},
        "subcutaneous": {"name": "Подкожная клетчатка (SQ)", "pO2": 30.0, "description": "Низкое давление O2. Доступно, но уязвимо для фиброза."},
        "extreme_hypoxia": {"name": "Зона фиброза / Выраженная гипоксия", "pO2": 10.0, "description": "Худший сценарий. Доступ кислорода сильно ограничен."}
    },
    "ibmir_kinetics": {
        "k_tf_thrombin": 0.08,
        "k_clot_platelet": 0.04,
        "heparin_inhibition": 0.15,
        "heparin_low_mw_inhibition": 0.08
    }
}

def load_parameters(filepath=None):
    """
    Loads parameter configuration from YAML file with fallback to built-in defaults.
    """
    target_path = filepath if filepath else DEFAULT_PARAMS_PATH
    if os.path.exists(target_path):
        try:
            with open(target_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                if isinstance(data, dict):
                    return data
        except Exception as e:
            print(f"[param_loader] Warning: Failed to parse {target_path} ({e}). Using fallbacks.")
    
    return DEFAULT_FALLBACK_PARAMS.copy()

def load_literature_parameters(filepath=None):
    """
    Loads full literature parameters database with citations.
    """
    target_path = filepath if filepath else LITERATURE_PARAMS_PATH
    if os.path.exists(target_path):
        try:
            with open(target_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                if isinstance(data, dict):
                    return data
        except Exception as e:
            print(f"[param_loader] Warning: Failed to parse {target_path} ({e}). Using empty literature dict.")
    return {}
