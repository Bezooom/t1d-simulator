import numpy as np
import trimesh
from skimage.measure import marching_cubes

class TPMSGenerator:
    """
    Класс для процедурной генерации трижды периодических минимальных поверхностей (TPMS)
    и расчета их биофизических параметров (пористость, площадь поверхности, SA/V).
    """
    def __init__(self, tpms_type="gyroid", box_size=1000.0, unit_cell=300.0, thickness=0.0, resolution=80):
        """
        tpms_type: тип минимальной поверхности ("gyroid" или "schwarz_p")
        box_size: размер кубического контейнера в микрометрах (мкм)
        unit_cell: размер элементарной ячейки в мкм (период)
        thickness: смещение изоповерхности t (управляет толщиной стенок гидрогеля)
        resolution: разрешение voxel-сетки (число точек по одной оси)
        """
        self.tpms_type = tpms_type.lower()
        self.box_size = box_size
        self.unit_cell = unit_cell
        self.thickness = thickness
        self.resolution = resolution
        
        self.mesh = None
        self.volume_fraction = 0.0
        self.porosity = 0.0
        self.surface_area = 0.0
        self.solid_volume = 0.0
        self.sav_ratio = 0.0
        
    def generate_voxel_grid(self):
        """
        Генерирует 3D сетку значений неявной функции TPMS.
        """
        x = np.linspace(0.0, self.box_size, self.resolution)
        y = np.linspace(0.0, self.box_size, self.resolution)
        z = np.linspace(0.0, self.box_size, self.resolution)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        
        # Перевод в фазовые координаты на основе элементарной ячейки
        kx = 2.0 * np.pi * X / self.unit_cell
        ky = 2.0 * np.pi * Y / self.unit_cell
        kz = 2.0 * np.pi * Z / self.unit_cell
        
        if self.tpms_type == "gyroid":
            # Неявное уравнение Гироида: sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x) - t
            F = np.sin(kx) * np.cos(ky) + np.sin(ky) * np.cos(kz) + np.sin(kz) * np.cos(kx) - self.thickness
        elif self.tpms_type == "schwarz_p":
            # Неявное уравнение Шварца P: cos(x) + cos(y) + cos(z) - t
            F = np.cos(kx) + np.cos(ky) + np.cos(kz) - self.thickness
        else:
            raise ValueError(f"Неизвестный тип TPMS: {self.tpms_type}. Поддерживаются 'gyroid' и 'schwarz_p'.")
            
        return F

    def build_mesh(self):
        """
        Экстрагирует полигональную сетку методом Marching Cubes с паддингом границ
        для создания замкнутого (водонепроницаемого) твердого тела.
        """
        F = self.generate_voxel_grid()
        
        # Трюк: оборачиваем voxel-массив отрицательными границами, 
        # чтобы marching cubes замкнул сетку на краях куба (watertight manifold).
        # Предполагается, что гидрогель находится в области F >= 0.0
        F_padded = np.pad(F, pad_width=1, mode='constant', constant_values=-10.0)
        
        # Шаг сетки в мкм
        spacing = self.box_size / (self.resolution - 1)
        
        try:
            verts, faces, normals, values = marching_cubes(
                F_padded, 
                level=0.0, 
                spacing=(spacing, spacing, spacing)
            )
            
            # Корректируем смещение вершин после паддинга (сдвиг на -1 шаг)
            verts -= spacing
            
            # Создаем trimesh объект
            self.mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals)
            
            # Дополнительное сглаживание и исправление нормалей (если они смотрят внутрь)
            trimesh.repair.fix_normals(self.mesh)
            trimesh.repair.fix_inversion(self.mesh)
            trimesh.repair.fix_winding(self.mesh)
            
            self.calculate_metrics()
            
        except Exception as e:
            # Если уровень изоповерхности вне диапазона значений сетки
            print(f"Ошибка экстракции 3D-сетки: {e}")
            self.mesh = None
            
        return self.mesh

    def calculate_metrics(self):
        """
        Рассчитывает геометрические параметры гидрогелевого каркаса.
        """
        if self.mesh is None:
            return
            
        # Общий объем кубического контейнера (мкм³)
        total_volume = self.box_size ** 3
        
        # Объем твердой фазы гидрогеля (мкм³)
        # В случае ошибок trimesh может вернуть отрицательный объем, берем модуль
        self.solid_volume = abs(self.mesh.volume)
        
        # Пористость (объемная доля пустот, доступная клеткам и жидкости)
        self.volume_fraction = self.solid_volume / total_volume
        self.porosity = (1.0 - self.volume_fraction) * 100.0
        
        # Площадь поверхности гидрогеля (мкм²)
        self.surface_area = self.mesh.area
        
        # Удельная площадь поверхности к объему твердой фазы (SA/V)
        # Переводим из мкм⁻¹ в см⁻¹ (1 мкм⁻¹ = 10^4 см⁻¹)
        if self.solid_volume > 0:
            self.sav_ratio = (self.surface_area / self.solid_volume) * 1.0e4
        else:
            self.sav_ratio = 0.0

    def export_stl(self, file_path):
        """
        Экспортирует сетку в файл STL.
        """
        if self.mesh is None:
            raise ValueError("Сетка не сгенерирована. Сначала вызовите build_mesh().")
            
        self.mesh.export(file_path, file_type="stl")
        print(f"Сетка успешно экспортирована в {file_path}")

    def get_stl_string(self):
        """
        Возвращает данные STL в виде ASCII строки.
        """
        if self.mesh is None:
            raise ValueError("Сетка не сгенерирована. Сначала вызовите build_mesh().")
            
        import mesh_generator
        return mesh_generator.export_to_stl_ascii(self.mesh.vertices, self.mesh.faces, solid_name=f"tpms_{self.tpms_type}")
