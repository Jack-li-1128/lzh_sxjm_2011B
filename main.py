import matplotlib.pyplot as plt
import random
import matplotlib.patheffects as path_effects
import matplotlib.cm as cm
import numpy as np
import matplotlib.colors as mcolors
import heapq
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster


def read_f_file(file_name):
    def line_count(data):
        line = 0
        while True:
            temp = data.readline()
            if temp == '':
                data.seek(0, 0)
                return line
            line += 1
    file = open(file_name, 'r', encoding='utf8')
    line = line_count(file)
    f = []
    for i in range(line):
        f.append(file.readline().strip('\n').split('\t'))
    return f


def main():
    f_point = read_f_file('data1.txt')
    f_relation = read_f_file('data2.txt')
    f_plant = read_f_file('data3.txt')
    f_inout_point = read_f_file('data4.txt')
    a_area_inout_point = [i[2] for i in f_inout_point]
    city_inout_point = [i[1] for i in f_inout_point]
    del a_area_inout_point[0]
    del city_inout_point[0]
    plant_name = [i[0] for i in f_plant]
    plant_point = [i[1] for i in f_plant]
    del plant_name[0]
    del plant_point[0]
    inout_point = a_area_inout_point + city_inout_point
    def distance(x1, y1, x2, y2):
        return ((x1-x2)**2 + (y1-y2)**2)**0.5

    def get_one_area_point(area_name=''):
        # points_list[[class_name,x,y],[class_name,x,y],......,[class_name,x,y]]
        points_dict = {}
        if area_name == '':
            for j in range(1, len(f_point)):
                points_dict[f_point[j][0]] = [float(f_point[j][1]), float(f_point[j][2]), f_point[j][3], 0]
                if f_point[j][0] in plant_point:
                    points_dict[f_point[j][0]].append(plant_name[plant_point.index(f_point[j][0])])
                else:
                    points_dict[f_point[j][0]].append('nor')
                if f_point[j][0] in inout_point:
                    points_dict[f_point[j][0]][3] = 1
                points_dict[f_point[j][0]].append(f_point[j][4])
        for j in range(1, len(f_point)):
            if f_point[j][3] == area_name:
                points_dict[f_point[j][0]] = [float(f_point[j][1]), float(f_point[j][2]), f_point[j][3], 0]
                if f_point[j][0] in plant_point:
                    points_dict[f_point[j][0]].append(plant_name[plant_point.index(f_point[j][0])])
                else:
                    points_dict[f_point[j][0]].append('nor')
                if f_point[j][0] in inout_point:
                    points_dict[f_point[j][0]][3] = 1
                points_dict[f_point[j][0]].append(f_point[j][4])
        return points_dict

    def create_road_point_relation(points_dict):
        # points_relation = {point_name1:[[distance11,point_name11],[distance11,point_name12],...],
        # point_name2:[[distance21,point_name11],[distance22,point_name12],...]...}
        points_relation = {}
        for i in range(1,len(f_relation)):
            try:
                points_dict.get(f_relation[i][1])
                points_dict.get(f_relation[i][0])
                try:
                    points_relation.get(f_relation[i][0])
                    temp_d = distance(points_dict[f_relation[i][0]][0], points_dict[f_relation[i][0]][1],
                                      points_dict[f_relation[i][1]][0], points_dict[f_relation[i][1]][1])
                    points_relation[f_relation[i][0]].append([temp_d, f_relation[i][1]])
                except:
                    temp_d = distance(points_dict[f_relation[i][0]][0], points_dict[f_relation[i][0]][1],
                                      points_dict[f_relation[i][1]][0], points_dict[f_relation[i][1]][1])
                    points_relation[f_relation[i][0]] = [[temp_d, f_relation[i][1]]]
            except:
                pass
        for i in range(1,len(f_relation)):
            try:
                points_dict.get(f_relation[i][1])
                points_dict.get(f_relation[i][0])
                try:
                    points_relation.get(f_relation[i][1])
                    temp_d = distance(points_dict[f_relation[i][0]][0], points_dict[f_relation[i][0]][1],
                                      points_dict[f_relation[i][1]][0], points_dict[f_relation[i][1]][1])
                    points_relation[f_relation[i][1]].append([temp_d, f_relation[i][0]])
                except:
                    temp_d = distance(points_dict[f_relation[i][0]][0], points_dict[f_relation[i][0]][1],
                                      points_dict[f_relation[i][1]][0], points_dict[f_relation[i][1]][1])
                    points_relation[f_relation[i][1]] = [[temp_d, f_relation[i][0]]]
            except:
                pass
        return points_relation

    def dijkstra(area, start, end):
        # 初始化最小堆
        graph = create_road_point_relation(get_one_area_point(area))
        queue = [(0, start, [])]
        seen = set()
        min_dist = {start: 0}
        while queue:
            (cost, v1, path) = heapq.heappop(queue)
            if v1 in seen:
                continue
            path = path + [v1]
            seen.add(v1)

            if v1 == end:
                return (cost, path)
            for c, v2 in graph.get(v1, []):
                if v2 in seen:
                    continue
                prev = min_dist.get(v2, None)
                next = cost + c
                if prev is None or next < prev:
                    min_dist[v2] = next
                    heapq.heappush(queue, (next, v2, path))
        return (float("inf"), [])

    def calculate_work(hold_area, area=''):
        nodes = get_one_area_point(area)
        w = {}
        for key, value in hold_area.items():
            t = 0
            for i in range(len(value)):
                for j, k in nodes.items():
                    if key in k:
                        pp = j
                        break
                dis, link = dijkstra(area, pp, value[i])
                t += float(nodes[value[i]][5]) * (dis * 0.2 + 30)
            w[key] = t
        return w

    def search(area=''):
        nodes = get_one_area_point(area)
        connections = create_road_point_relation(nodes)
        plant_reach_area = dict()
        max_reach = 30
        def track_back_search_area(point, dis, visited):
            for j in connections.get(point, []):
                new_dis = dis + j[0]
                if j[1] not in visited and new_dis < max_reach:
                    visited.add(j[1])
                    plant_reach.add(j[1])
                    track_back_search_area(j[1], new_dis, visited)
                    visited.remove(j[1])  # 回溯时移除节点
        for i in range(len(plant_point)):
            if connections.get(plant_point[i]) != None:
                plant_reach = set([plant_point[i]])
                dis = 0
                track_back_search_area(plant_point[i], dis, set([plant_point[i]]))
                plant_reach_area[plant_name[i]] = plant_reach
        return plant_reach_area

    def search_by_step(area):
        nodes = get_one_area_point(area)
        connections = create_road_point_relation(nodes)
        plant_hold = {}
        all_point = []
        for key, value in nodes.items():
            if value[4] != 'nor':
                plant_hold[value[4]] = []
            all_point.append(key)
        all_point = set(all_point)
        rate = [1 for i in range(len(plant_hold))]
        max_reach = [(1 / 2) * rate[j] ** 2 for j in range(len(plant_hold))]
        for i in range(5, 240):
            plant_reach_area = dict()
            max_reach = [(1/2) * rate[j] ** 2 + max_reach[j] for j in range(len(plant_hold))]
            def track_back_search_area(point, dis, visited):
                for k in connections.get(point, []):
                    new_dis = dis + k[0]
                    if k[1] not in visited and new_dis < max_reach[j]:
                        visited.add(k[1])
                        plant_reach.add(k[1])
                        track_back_search_area(k[1], new_dis, visited)
                        visited.remove(k[1])  # 回溯时移除节点
            for j in range(len(plant_point)):
                if connections.get(plant_point[j]) != None:
                    plant_reach = set([plant_point[j]])
                    dis = 0
                    if max_reach[j] < 35:
                        track_back_search_area(plant_point[j], dis, set([plant_point[j]]))
                        plant_reach_area[plant_name[j]] = plant_reach
            for key,value in plant_reach_area.items():
                for p in value:
                    if p in all_point:
                        all_point -= {p}
                        plant_hold[key].append(p)
            c_t = calculate_work(plant_hold, area)
            av = [value for k, value in c_t.items()]
            rate = [(sum(av)/len(c_t))/value for k, value in c_t.items()]
        print(all_point)
        print(len(all_point))
        print(plant_hold)
        print(calculate_work(plant_hold))
        return plant_hold


    def draw_map(area=''):
        color_map = {'A': 'red', 'B': 'blue', 'C': 'green', 'D': 'yellow', 'E': 'purple'}  # 可以根据需要添加更多类别和颜色
        fig, ax = plt.subplots(figsize=(50, 50), dpi=300)
        nodes = get_one_area_point(area)
        # 绘制节点
        for node, (x, y, category, entrance, description, p) in nodes.items():
            color = color_map.get(category, 'black')  # 如果类别不在映射中，默认为黑色
            ax.plot(x, y, 'o', color=color)
            ax.text(x, y, node, fontsize=12, ha='right')
            #label = f'*{description}' if entrance == 1 else description
            #ax.text(x, y, label, fontsize=10, ha='right')
        connections = create_road_point_relation(nodes)
        # 绘制连接关系
        for node, edges in connections.items():
            x1, y1, _, _, _, _ = nodes[node]
            for edge in edges:
                distance, target = edge
                if target in nodes:
                    x2, y2, _, _, _, _ = nodes[target]
                    ax.plot([x1, x2], [y1, y2], 'k-', lw=0.5)  # 节点之间的线
                    mid_x = (x1 + x2) / 2
                    mid_y = (y1 + y2) / 2
                    # 在线的中点处标注距离
                    ax.text(mid_x, mid_y, f'{distance / 10:.1f}', fontsize=8, ha='center', va='center', color='blue')
        collections = search_by_step(area)
        all_point_relate_plant = {}
        collection_colors = colors = [
            "#FF0000",  # 红色
            "#00FF00",  # 绿色
            "#0000FF",  # 蓝色
            "#FFFF00",  # 黄色
            "#FF00FF",  # 品红
            "#00FFFF",  # 青色
            "#FFA500",  # 橙色
            "#800080",  # 紫色
            "#00FF7F",  # 春绿
            "#FF1493",  # 深粉红
            "#1E90FF",  # 道奇蓝
            "#32CD32",  # 石灰绿
            "#FF4500",  # 橙红
            "#8A2BE2",  # 蓝紫色
            "#ADFF2F",  # 绿黄色
            "#FF69B4",  # 热粉红
            "#40E0D0",  # 绿松石
            "#FF6347",  # 番茄
            "#EE82EE",  # 紫罗兰
            "#7FFF00",   # 查特酒绿
            "#FF0000",  # 红色
            "#00FF00",  # 绿色
            "#0000FF",  # 蓝色
            "#FF0000",  # 红色
            "#00FF00",  # 绿色
            "#0000FF",  # 蓝色
            "#FFFF00",  # 黄色
            "#FF00FF",  # 品红
            "#00FFFF",  # 青色
            "#FFA500",  # 橙色
            "#800080",  # 紫色
            "#00FF7F",  # 春绿
            "#FF1493",  # 深粉红
            "#1E90FF",  # 道奇蓝
            "#32CD32",  # 石灰绿
            "#FF4500",  # 橙红
            "#8A2BE2",  # 蓝紫色
            "#ADFF2F",  # 绿黄色
            "#FF69B4",  # 热粉红
            "#40E0D0",  # 绿松石
            "#FF6347",  # 番茄
            "#EE82EE",  # 紫罗兰
            "#7FFF00",  # 查特酒绿
            "#FF0000",  # 红色
            "#00FF00",  # 绿色
            "#0000FF",  # 蓝色
        ]
        for key, value in collections.items():
            for i in value:
                try:
                    all_point_relate_plant[i].append(key)
                except:
                    all_point_relate_plant[i] = [key]
        for point_name, collection_plants in all_point_relate_plant.items():
            x, y, _, _, class_n, _ = nodes[point_name]
            if class_n != 'nor':
                ax.text(x, y, class_n, fontsize=12, ha='left', va='center',
                            color=collection_colors[(int(class_n[1:]) - 1)%20])
            for i in range(len(collection_plants)):
                ax.scatter(x, y, s=(i + 1) * 100, edgecolor=collection_colors[(int(collection_plants[i][1:]) - 1)%20],
                               facecolor='none', linewidth=0.5)
        # 设置图形范围
        x_values, y_values, _, _, _, _ = zip(*nodes.values())
        ax.set_xlim(min(x_values) - 10, max(x_values) + 10)
        ax.set_ylim(min(y_values) - 10, max(y_values) + 10)
        # 添加标签和标题
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Map')
        ax.legend()
        # 保存图形
        plt.savefig('map.png')

    #print(search_by_step('A'))
    #draw_map('')

    def found_closed_dis(area):
        node = get_one_area_point(area)
        ptp = {}
        for key, value in node.items():
            if value[3] == 1:
                ptp[key] = []
                for i in range(20):
                    t = dijkstra(area, key, str(i+1))
                    ptp[key].append(['A'+str(i+1),t[0]])
        return ptp

    def find_min_sum_distance(graph):
        # 提取点的名称和距离
        points = list(graph.keys())
        all_connections = []
        for point in points:
            all_connections.extend([conn[0] for conn in graph[point]])
        all_connections = list(set(all_connections))
        # 构建成本矩阵
        cost_matrix = np.zeros((len(points), len(all_connections)))
        for i, point in enumerate(points):
            for j, conn in enumerate(all_connections):
                for connection in graph[point]:
                    if connection[0] == conn:
                        cost_matrix[i][j] = connection[1]
                        break
                else:
                    cost_matrix[i][j] = np.inf  # 如果没有连接，设置为无穷大
        # 使用匈牙利算法求解最小成本匹配
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        # 计算最小距离和
        min_distance_sum = cost_matrix[row_ind, col_ind].sum()
        # 获取匹配结果
        matching = [(points[i], all_connections[j]) for i, j in zip(row_ind, col_ind)]
        return min_distance_sum, matching
    '''
    min_distance_sum, matching = find_min_sum_distance(found_closed_dis('A'))
    print(f"最小距离和是: {min_distance_sum}")
    print("匹配结果是:")
    for match in matching:
        print(match)
    '''

    def release_plant(area=''):
        nodes = get_one_area_point(area)
        connections = create_road_point_relation(nodes)
        point_reach_area = dict()
        max_reach = 30
        def track_back_search_area(point, dis, visited):
            for j in connections.get(point, []):
                new_dis = dis + j[0]
                if j[1] not in visited and new_dis < max_reach:
                    visited.add(j[1])
                    plant_reach.add(j[1])
                    track_back_search_area(j[1], new_dis, visited)
                    visited.remove(j[1])  # 回溯时移除节点
        for key,value in nodes.items():
            plant_reach = set([key])
            dis = 0
            track_back_search_area(key, dis, set(key))
            point_reach_area[key] = list(plant_reach)
        work = calculate_work(point_reach_area, '')
        plant = []
        for i in range(len(work)):
            print(work.values())
            if len(work) != 0 and len(plant) < 80:
                #max_key = max(work, key=work.get)
                average_value = (sum(work.values()) / len(work))*2
                # 找到离平均值最近的value对应的key
                closest_key = min(work, key=lambda k: abs(work[k] - average_value))
                plant.append(closest_key)
                for j in point_reach_area[closest_key]:
                    if j in work.keys():
                        del work[j]
            else:
                break
        return plant
    #print(release_plant(''))

    def found_crime(area):
        nodes = get_one_area_point(area)
        connections = create_road_point_relation(nodes)
        plant_reach_area = dict()
        def track_back_search_area(point, dis, visited):
            for j in connections.get(point, []):
                new_dis = dis + j[0]
                if j[1] not in visited and new_dis < 30:
                    visited.add(j[1])
                    plant_reach.add(j[1])
                    track_back_search_area(j[1], new_dis, visited)
                    visited.remove(j[1])  # 回溯时移除节点
        for i in range(30,32):
            plant_reach = set(['32'])
            dis = 0
            track_back_search_area('32', dis, set(['32']))
            plant_reach_area['32'] = plant_reach
            print(plant_reach_area)
            all_connected_points = plant_reach_area['32']
            round_p = set(all_connected_points)
            for j in list(all_connected_points):
                for k in connections.get(j, []):
                    round_p.add(k[1])
            print(round_p)
            round_p -= all_connected_points
            print(round_p)
        return






main()

