import numpy as np

def generate_box_mesh(L_microns, width_microns=5000, height_microns=10000):
    """
    Генерирует 3D полигональную сетку плоской пластины (Box).
    Размеры: 2L x width x height
    """
    x = L_microns
    y = width_microns / 2.0
    z = height_microns / 2.0
    
    # 8 вершин коробки
    vertices = np.array([
        [-x, -y, -z],
        [ x, -y, -z],
        [ x,  y, -z],
        [-x,  y, -z],
        [-x, -y,  z],
        [ x, -y,  z],
        [ x,  y,  z],
        [-x,  y,  z]
    ])
    
    # 12 треугольных граней (индексы вершин)
    faces = np.array([
        # Низ
        [0, 2, 1], [0, 3, 2],
        # Верх
        [4, 5, 6], [4, 6, 7],
        # Спереди
        [0, 1, 5], [0, 5, 4],
        # Справа
        [1, 2, 6], [1, 6, 5],
        # Сзади
        [2, 3, 7], [2, 7, 6],
        # Слева
        [3, 0, 4], [3, 4, 7]
    ])
    
    return vertices, faces

def generate_cylinder_mesh(R_microns, length_microns=10000, num_segments=32):
    """
    Генерирует 3D сетку цилиндра.
    Радиус: R, Длина: length
    """
    vertices = []
    faces = []
    
    half_len = length_microns / 2.0
    
    # 1. Генерируем вершины нижнего и верхнего оснований
    for z in [-half_len, half_len]:
        for i in range(num_segments):
            theta = 2.0 * np.pi * i / num_segments
            x = R_microns * np.cos(theta)
            y = R_microns * np.sin(theta)
            vertices.append([x, y, z])
            
    # Добавляем центральные точки оснований для закрытия крышек
    center_bottom_idx = len(vertices)
    vertices.append([0.0, 0.0, -half_len])
    
    center_top_idx = len(vertices)
    vertices.append([0.0, 0.0, half_len])
    
    vertices = np.array(vertices)
    
    # 2. Боковая поверхность
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        
        # Индексы вершин на нижнем и верхнем круге
        b1 = i
        b2 = next_i
        t1 = i + num_segments
        t2 = next_i + num_segments
        
        # Два треугольника на каждый прямоугольный сегмент стенки
        faces.append([b1, b2, t2])
        faces.append([b1, t2, t1])
        
    # 3. Нижняя крышка (кап)
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        faces.append([next_i, i, center_bottom_idx])
        
    # 4. Верхняя крышка (кап)
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        faces.append([i, next_i, center_top_idx])
        
    return vertices, np.array(faces)

def generate_sphere_mesh(R_microns, num_latitude=16, num_longitude=32):
    """
    Генерирует 3D сетку сферы (UV Sphere).
    """
    vertices = []
    faces = []
    
    # Вершины
    for i in range(num_latitude + 1):
        theta = np.pi * i / num_latitude
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        
        for j in range(num_longitude):
            phi = 2.0 * np.pi * j / num_longitude
            x = R_microns * sin_theta * np.cos(phi)
            y = R_microns * sin_theta * np.sin(phi)
            z = R_microns * cos_theta
            vertices.append([x, y, z])
            
    vertices = np.array(vertices)
    
    # Грани
    for i in range(num_latitude):
        for j in range(num_longitude):
            next_j = (j + 1) % num_longitude
            
            # Индексы вершин сетки
            # Строка i, колонка j
            p1 = i * num_longitude + j
            p2 = i * num_longitude + next_j
            p3 = (i + 1) * num_longitude + j
            p4 = (i + 1) * num_longitude + next_j
            
            # Строим треугольники (исключая вырожденные на полюсах)
            if i > 0:
                faces.append([p1, p2, p4])
            if i < num_latitude - 1:
                faces.append([p1, p4, p3])
                
    return vertices, np.array(faces)

def export_to_stl_ascii(vertices, faces, solid_name="encapsulated_graft"):
    """
    Преобразует вершины и грани в текстовый (ASCII) формат STL.
    """
    lines = [f"solid {solid_name}"]
    
    for face in faces:
        v1 = vertices[face[0]]
        v2 = vertices[face[1]]
        v3 = vertices[face[2]]
        
        # Векторное произведение для расчета нормали
        val1 = v2 - v1
        val2 = v3 - v1
        normal = np.cross(val1, val2)
        norm = np.linalg.norm(normal)
        if norm > 0:
            normal = normal / norm
        else:
            normal = np.array([0.0, 0.0, 0.0])
            
        lines.append(f"  facet normal {normal[0]:.6e} {normal[1]:.6e} {normal[2]:.6e}")
        lines.append("    outer loop")
        lines.append(f"      vertex {v1[0]:.6f} {v1[1]:.6f} {v1[2]:.6f}")
        lines.append(f"      vertex {v2[0]:.6f} {v2[1]:.6f} {v2[2]:.6f}")
        lines.append(f"      vertex {v3[0]:.6f} {v3[1]:.6f} {v3[2]:.6f}")
        lines.append("    endloop")
        lines.append("  endfacet")
        
    lines.append(f"endsolid {solid_name}")
    return "\n".join(lines)
