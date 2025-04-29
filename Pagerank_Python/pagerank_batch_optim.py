length = 9500

# 我们利用及时计算的思想来处理：
hash_dict = {}
is_calculated2 = [False] * (length)


def HashMap(node_id):
    return hash_dict[node_id]

def ReverseHashMap(hash_value):
    for key, value in hash_dict.items():
        if value == hash_value:
            return key

def converge(r_old, r_new, epsilon):
    sum = 0
    for i in range(len(r_old)):
        sum += abs(r_new[i] - r_old[i])
    #print("sum:", sum)
    return sum >= epsilon  # 返回True表示继续迭代

    

def pagerank_all(epsilon=0.0001, beta=0.85, file_path='Data.txt', total_batch=15):
    global is_calculated2
    r_old = [1/length] * length
    round = 0
    
    while True:
        round += 1
        print("round:", round)
        r_new = [(1 - beta) / length] * length
        is_calculated2 = [False] * length

        for batch_id in range(total_batch):  # 每轮都读取所有 batch！
            line_num = 0
            total_lines = 150000
            batch_size = total_lines // total_batch
            transM_instant = {}
            with open(file_path, "r") as file:
                for line in file:
                    line_num += 1
                    if batch_id*batch_size <= line_num and line_num < (batch_id+1)*batch_size:
                        parts = line.strip().split()
                        from_node_id = int(parts[0])
                        to_node_id = int(parts[1])
                        if is_calculated2[HashMap(from_node_id)] == False and from_node_id not in transM_instant: #需要hash
                            transM_instant[from_node_id] = {"degree": 1, "dest": [to_node_id]}
                        elif is_calculated2[HashMap(from_node_id)] == False and from_node_id in transM_instant:
                            transM_instant[from_node_id]["degree"] += 1
                            transM_instant[from_node_id]["dest"].append(to_node_id)

                        if is_calculated2[HashMap(to_node_id)] == False and to_node_id not in transM_instant: #需要hash
                            transM_instant[to_node_id] = {"degree": 0, "dest": []}
                        
                    else:
                        continue

            line_num = 0
            with open(file_path, "r") as file:
                for line2 in file:
                    line_num += 1
                    if batch_id*batch_size  <= line_num and line_num < (batch_id+1)*batch_size:
                        continue
                    parts2 = line2.strip().split()
                    from_node_id2 = int(parts2[0])
                    to_node_id2 = int(parts2[1])
                    if from_node_id2 in transM_instant:
                        transM_instant[from_node_id2]["degree"] += 1
                        transM_instant[from_node_id2]["dest"].append(to_node_id2)
            
            if transM_instant == {}:
                continue
            for key in transM_instant:
                is_calculated2[HashMap(key)] = True
            #print("batch:", batch_id + 1)
            for source in transM_instant:
                degree = transM_instant[source]["degree"]
                if degree == 0:
                    contribution = (beta) * r_old[HashMap(source)] / length
                    for k in range(length):
                        r_new[k] += contribution
                else:
                    contribution = beta * r_old[HashMap(source)] / degree
                    for dest in transM_instant[source]["dest"]:
                        j = HashMap(dest)
                        r_new[j] += contribution

        if not converge(r_old, r_new, epsilon):
            break
        r_old = r_new
    
    return r_new

        
def GetTotalHash(file_path):
    global hash_dict
    with open(file_path, "r") as file:
        dict_num = 0
        for line in file:
            parts = line.strip().split()
            from_node_id = int(parts[0])
            to_node_id = int(parts[1])
            if from_node_id not in hash_dict:
                hash_dict[from_node_id] = dict_num
                dict_num += 1
            if to_node_id not in hash_dict:
                hash_dict[to_node_id] = dict_num
                dict_num += 1
                
GetTotalHash("Data.txt")
vector_r = pagerank_all(total_batch=150)

def get_top_100_indices(r_new):
    indexed_values = [(index, value) for index, value in enumerate(r_new)]
    indexed_values.sort(key=lambda x: x[1], reverse=True)
    top_100_indices = [
        (ReverseHashMap(index), value) for index, value in indexed_values[:100]
    ]
    return top_100_indices


top_100 = get_top_100_indices(vector_r)
print(top_100)