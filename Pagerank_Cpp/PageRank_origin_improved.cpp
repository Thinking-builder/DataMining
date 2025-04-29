#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iomanip>
using namespace std;

const int MAX_LENGTH = 9500;       // 点数上限
const double BETA = 0.85;
const double EPSILON = 0.0001;
const int TOTAL_BATCH = 1;         // 批次数
const int BATCH_LINES = 150000 / TOTAL_BATCH; // 每批次行数

// 哈希映射结构
int hash_map[MAX_LENGTH * 10];         // node_id -> index
int reverse_map[MAX_LENGTH];          // index -> node_id
bool used_id[MAX_LENGTH * 10];
int current_index = 0;

// 辅助结构
int degree[MAX_LENGTH];               // 节点出度
double r_old[MAX_LENGTH];             // 上一轮 PageRank 分数
double r_new[MAX_LENGTH];             // 当前轮 PageRank 分数

int get_hash(int node_id) {
    if (!used_id[node_id]) {
        used_id[node_id] = true;
        hash_map[node_id] = current_index;
        reverse_map[current_index] = node_id;
        return current_index++;
    }
    return hash_map[node_id];
}

// 检查收敛
bool converge(double* r_old, double* r_new) {
    double sum = 0.0;
    for (int i = 0; i < current_index; ++i) {
        sum += fabs(r_new[i] - r_old[i]);
    }
    return sum < EPSILON;
}

// 获取所有节点并计算出度
void GetTotalHashAndDegree(const char* file_path) {
    ifstream fin(file_path);
    int from_id, to_id;
    while (fin >> from_id >> to_id) {
        int from_idx = get_hash(from_id);
        int to_idx = get_hash(to_id);
        degree[from_idx]++; // 计算出度
    }
    fin.close();
}

// PageRank 计算
void pagerank_all(const char* file_path) {
    // 初始化 PageRank 分数
    for (int i = 0; i < current_index; ++i) {
        r_old[i] = 1.0 / current_index;
    }
    bool cont = true;
    while (cont) {
        // 初始化 r_new
        for (int i = 0; i < current_index; ++i) {
            r_new[i] = (1.0 - BETA) / current_index;
        }
        // 读取文件并按批次处理
        ifstream fin(file_path);
        char line[256];
        int batch_line_count = 0;
        int batch_id = 0;

        while (fin.getline(line, sizeof(line))) {
            int from, to;
            if (sscanf(line, "%d %d", &from, &to) == 2) {
                int from_idx = get_hash(from);
                int to_idx = get_hash(to);
                // 添加贡献
                if (degree[from_idx] > 0) {
                    r_new[to_idx] += BETA * (r_old[from_idx] / degree[from_idx]);
                }

                // 分批次处理
                if (++batch_line_count >= BATCH_LINES && !fin.eof()) {
                    batch_line_count = 0;
                    if (++batch_id >= TOTAL_BATCH){
                        break;
                    }
                }
            }
        }
        fin.close();
        // 处理悬挂节点（dangling nodes）
        double M = 0.0;
        for (int i = 0; i < current_index; ++i) {
            if (degree[i] == 0) {
                M += r_old[i];
            }
        }
        for (int i = 0; i < current_index; ++i) {
            r_new[i] += BETA * (M / current_index);
        }

        // 检查收敛
        cont = !converge(r_old, r_new);
        // 更新 r_old
        for (int i = 0; i < current_index; ++i) {
            r_old[i] = r_new[i];
        }
    }

    // 输出前 100 个节点
    ofstream ofs("top100.txt");
    if (!ofs.is_open()) {
        cerr << "无法打开输出文件。" << endl;
        return;
    }

    for (int i = 0; i < 100 && i < current_index; ++i) {
        int max_idx = 0;
        for (int j = 1; j < current_index; ++j) {
            if (r_old[j] > r_old[max_idx]) max_idx = j;
        }
        printf("Node %d: %.8lf\n", reverse_map[max_idx], r_old[max_idx]);
        ofs << "Node " << reverse_map[max_idx] << ": " << fixed << setprecision(8) << r_old[max_idx] << endl;
        r_old[max_idx] = -1; // 防止重复
    }
    ofs.close();
}

int main() {
    memset(used_id, 0, sizeof(used_id));
    memset(degree, 0, sizeof(degree));
    GetTotalHashAndDegree("Data.txt");
    pagerank_all("Data.txt");
    return 0;
}