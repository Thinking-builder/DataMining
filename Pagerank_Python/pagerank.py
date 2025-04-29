transM = {}


def read_data(file_path):
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split()
            from_node_id = int(parts[0])
            to_node_id = int(parts[1])
            if from_node_id not in transM:
                transM[from_node_id] = {"degree": 1, "dest": [to_node_id]}
            else:
                transM[from_node_id]["degree"] += 1
                transM[from_node_id]["dest"].append(to_node_id)
            if to_node_id not in transM:
                transM[to_node_id] = {"degree": 0, "dest": []}


file_path = "Data.txt"
read_data(file_path)

import numpy as np

# from time import sleep

# sleep(10)


def converge(r_old, r_new, epsilon):
    return np.sum(np.abs(r_new - r_old)) >= epsilon  # 返回True表示继续迭代


HashDict = {}
# 假设transM已定义且包含所有节点信息
ii = 0
for key in transM:
    HashDict[key] = ii
    ii += 1
length = len(HashDict)  # 确保length正确初始化


def HashMap(Matrixkey):
    return HashDict[Matrixkey]


def ReverseHashMap(index):
    for key, value in HashDict.items():
        if value == index:
            return key


def pagerank(Matrix, epsilon=0.0001, beta=0.85):
    r_old = np.array([1.0 / length] * length)
    iteration_count = 0  # 可选：跟踪迭代次数

    while True:
        r_new = np.array(
            [(1 - beta) / length] * length, dtype=np.float64
        )  # 每次迭代重置基值

        for i in range(length):
            source = ReverseHashMap(i)
            degree = Matrix[source]["degree"]

            if degree == 0:
                # 处理悬挂节点：贡献均分给所有节点
                contribution = beta * r_old[i] / length
                r_new += contribution
            else:
                contribution = beta * r_old[i] / degree
                for dest in Matrix[source]["dest"]:
                    j = HashMap(dest)
                    r_new[j] += contribution

        if not converge(r_old, r_new, epsilon):
            break

        r_old = r_new.copy()
        iteration_count += 1

    return r_new


def get_top_100_indices(r_new):
    indexed_values = [(index, value) for index, value in enumerate(r_new)]
    indexed_values.sort(key=lambda x: x[1], reverse=True)
    top_100_indices = [
        (ReverseHashMap(index), value) for index, value in indexed_values[:100]
    ]
    return top_100_indices


top_100 = get_top_100_indices(pagerank(transM))
print(top_100)
