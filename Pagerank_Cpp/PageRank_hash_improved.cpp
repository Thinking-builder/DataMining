// 之前的一个版本，使用了线性探测法的哈希表，
// #include <cstdio>
// #include <cstdlib>
// #include <cstdint>
// #include <cmath>
// #include <cstring>

// const double BETA = 0.85;
// const double EPSILON = 1e-4;

// // 哈希表元素
// struct HashNode {
//     int key; // node_id
//     int val; // index
// };

// int hash_size;
// HashNode *hash_table;
// int N = 0; // 实际节点数
// int *idx2id;
// uint16_t *degree;
// float *r_old, *r_new;

// // 简单哈希函数
// inline int hash_func(int key) {
//     return (key * 23333u) % hash_size;
// }

// // 哈希插入 key -> N
// void insert(int key) {
//     int h = hash_func(key);
//     while (hash_table[h].key != -1) {
//         if (hash_table[h].key == key) return; // 已存在
//         h = (h + 1) % hash_size;
//     }
//     hash_table[h].key = key;
//     hash_table[h].val = N++;
// }

// // 哈希查询 key -> index
// int find(int key) {
//     int h = hash_func(key);
//     while (hash_table[h].key != -1) {
//         if (hash_table[h].key == key) return hash_table[h].val;
//         h = (h + 1) % hash_size;
//     }
//     return -1;
// }

// int main() {
//     FILE *fin = fopen("Data.txt", "r");
//     if (!fin) return 1;

//     // 预估节点数，开大一点防止冲突
//     hash_size = 40000; // 2倍左右  (根据实际数据可以调整)
//     hash_table = (HashNode*)malloc(sizeof(HashNode) * hash_size);
//     for (int i = 0; i < hash_size; ++i) hash_table[i].key = -1;

//     int u, v;
//     // 第一遍：建立节点哈希表
//     while (fscanf(fin, "%d %d", &u, &v) == 2) {
//         insert(u);
//         insert(v);
//     }

//     // 分配空间
//     idx2id = (int*)malloc(N * sizeof(int));
//     degree = (uint16_t*)calloc(N, sizeof(uint16_t));
//     r_old = (float*)malloc(N * sizeof(float));
//     r_new = (float*)malloc(N * sizeof(float));

//     if (!idx2id || !degree || !r_old || !r_new) return 1;

//     // 第二遍：记录idx2id
//     for (int i = 0; i < hash_size; ++i) {
//         if (hash_table[i].key != -1) {
//             idx2id[hash_table[i].val] = hash_table[i].key;
//         }
//     }

//     // 第三遍：统计出度
//     rewind(fin);
//     while (fscanf(fin, "%d %d", &u, &v) == 2) {
//         int ui = find(u);
//         if (ui >= 0) degree[ui]++;
//     }

//     // 初始化PageRank
//     for (int i = 0; i < N; ++i) r_old[i] = 1.0f / N;

//     // PageRank迭代
//     while (1) {
//         for (int i = 0; i < N; ++i) r_new[i] = (1.0f - BETA) / N;
//         float dsum = 0;
//         for (int i = 0; i < N; ++i) if (degree[i] == 0) dsum += r_old[i];
//         rewind(fin);
//         while (fscanf(fin, "%d %d", &u, &v) == 2) {
//             int ui = find(u), vi = find(v);
//             if (ui >= 0 && vi >= 0 && degree[ui]) {
//                 r_new[vi] += BETA * (r_old[ui] / degree[ui]);
//             }
//         }
//         float avg = BETA * dsum / N;
//         for (int i = 0; i < N; ++i) r_new[i] += avg;
//         float diff = 0;
//         for (int i = 0; i < N; ++i) diff += fabsf(r_new[i] - r_old[i]);
//         if (diff < EPSILON) break;
//         memcpy(r_old, r_new, N * sizeof(float));
//     }
//     fclose(fin);

//     // Top-100 插入排序
//     int K = N < 100 ? N : 100;
//     int *top_id = (int*)malloc(K * sizeof(int));
//     float *top_val = (float*)malloc(K * sizeof(float));
//     int cnt = 0;
//     for (int i = 0; i < N; ++i) {
//         float v = r_old[i];
//         int k = cnt;
//         if (cnt < K) {
//             top_id[cnt] = i; top_val[cnt++] = v;
//         } else if (v > top_val[0]) {
//             top_val[0] = v; top_id[0] = i; k = 0;
//         } else continue;
//         while (k + 1 < cnt && top_val[k] > top_val[k+1]) {
//             float tv = top_val[k+1]; top_val[k+1] = top_val[k]; top_val[k] = tv;
//             int ti = top_id[k+1]; top_id[k+1] = top_id[k]; top_id[k] = ti;
//             ++k;
//         }
//     }

//     // 输出
//     FILE *ofs = fopen("top100.txt", "w");
//     for (int i = cnt - 1; i >= 0; --i) {
//         int oid = idx2id[top_id[i]];
//         fprintf(ofs, "Node %d: %.8f\n", oid, top_val[i]);
//         printf("Node %d: %.8f\n", oid, top_val[i]);
//     }
//     fclose(ofs);

//     // 释放内存
//     free(hash_table);
//     free(idx2id);
//     free(degree);
//     free(r_old);
//     free(r_new);
//     free(top_id);
//     free(top_val);

//     return 0;
// }

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>

static const double BETA = 0.85;
static const double EPS   = 1e-4;

// 哈希表大小，2^15 = 32768，需大于 2*N
static const int HT_BITS = 15;
static const int HT_SIZE = 1 << HT_BITS;
static const uint32_t EMPTY = 0xFFFFFFFFu;

// 哈希映射函数
static inline uint32_t hash32(uint32_t x) {
    // Knuth's multiplicative hash
    return (x * 2654435761u) >> (32 - HT_BITS);
}

int main() {
    FILE* fin = fopen("Data.txt", "r");
    if (!fin) return 1;

    // 1. 构建哈希表收集唯一节点ID
    uint32_t* ht_keys = (uint32_t*)malloc(sizeof(uint32_t) * HT_SIZE);
    uint16_t* ht_vals = (uint16_t*)malloc(sizeof(uint16_t) * HT_SIZE);
    for (int i = 0; i < HT_SIZE; ++i) ht_keys[i] = EMPTY;

    uint32_t u, v;
    int N = 0;
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        // 插入 u
        uint32_t h = hash32(u);
        while (ht_keys[h] != EMPTY && ht_keys[h] != u) h = (h + 1) & (HT_SIZE - 1);
        if (ht_keys[h] == EMPTY) { ht_keys[h] = u; ++N; }
        // 插入 v
        h = hash32(v);
        while (ht_keys[h] != EMPTY && ht_keys[h] != v) h = (h + 1) & (HT_SIZE - 1);
        if (ht_keys[h] == EMPTY) { ht_keys[h] = v; ++N; }
    }

    // 2. 提取、排序、索引映射
    uint32_t* ids = (uint32_t*)malloc(sizeof(uint32_t) * N);
    int idx = 0;
    for (int i = 0; i < HT_SIZE; ++i) {
        if (ht_keys[i] != EMPTY) ids[idx++] = ht_keys[i];
    }
    // 排序
    qsort(ids, N, sizeof(uint32_t), [](const void* a, const void* b){
        uint32_t A = *(const uint32_t*)a;
        uint32_t B = *(const uint32_t*)b;
        return (A > B) - (A < B);
    });
    // 重构哈希表值：id -> 索引
    for (int i = 0; i < HT_SIZE; ++i) ht_keys[i] = EMPTY;
    for (int i = 0; i < N; ++i) {
        uint32_t key = ids[i];
        uint32_t h = hash32(key);
        while (ht_keys[h] != EMPTY && ht_keys[h] != key) h = (h + 1) & (HT_SIZE - 1);
        ht_keys[h] = key;
        ht_vals[h] = (uint16_t)i;
    }

    // 3. 统计出度
    uint16_t* degree = (uint16_t*)calloc(N, sizeof(uint16_t));
    rewind(fin);
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        uint32_t h = hash32(u);
        while (ht_keys[h] != u) h = (h + 1) & (HT_SIZE - 1);
        degree[ ht_vals[h] ]++;
    }

    // 4. 分配与初始化 PageRank 数组
    float* r_old = (float*)malloc(sizeof(float) * N);
    float* r_new = (float*)malloc(sizeof(float) * N);
    for (int i = 0; i < N; ++i) r_old[i] = 1.0f / N;

    // 5. PageRank 迭代（重读文件）
    while (true) {
        double dsum = 0.0;
        for (int i = 0; i < N; ++i) {
            if (degree[i] == 0) dsum += r_old[i];
            r_new[i] = (1.0f - BETA) / N;
        }
        rewind(fin);
        while (fscanf(fin, "%u %u", &u, &v) == 2) {
            uint32_t h = hash32(u);
            while (ht_keys[h] != u) h = (h + 1) & (HT_SIZE - 1);
            int ui = ht_vals[h];
            if (!degree[ui]) continue;
            h = hash32(v);
            while (ht_keys[h] != v) h = (h + 1) & (HT_SIZE - 1);
            int vi = ht_vals[h];
            r_new[vi] += (float)(BETA * (r_old[ui] / degree[ui]));
        }
        float avg = (float)(BETA * dsum / N);
        double diff = 0.0;
        for (int i = 0; i < N; ++i) {
            r_new[i] += avg;
            diff += fabsf(r_new[i] - r_old[i]);
        }
        if (diff < EPS) break;
        float* tmp = r_old; r_old = r_new; r_new = tmp;
    }
    fclose(fin);
    free(degree);
    free(r_new);

    // 6. 在线 Top-100
    int K = N < 100 ? N : 100;
    float  top_val[100];
    uint32_t top_id[100];
    int cnt_top = 0;
    for (int i = 0; i < N; ++i) {
        float val = r_old[i];
        uint32_t id = ids[i];
        if (cnt_top < K) {
            int j = cnt_top++;
            top_val[j] = val; top_id[j] = id;
            while (j > 0 && top_val[j-1] > top_val[j]) {
                // 交换
                float tv = top_val[j]; top_val[j] = top_val[j-1]; top_val[j-1] = tv;
                uint32_t ti = top_id[j]; top_id[j] = top_id[j-1]; top_id[j-1] = ti;
                --j;
            }
        } else if (val > top_val[0]) {
            top_val[0] = val; top_id[0] = id;
            int j = 0;
            while (j + 1 < K && top_val[j] > top_val[j+1]) {
                float tv = top_val[j]; top_val[j] = top_val[j+1]; top_val[j+1] = tv;
                uint32_t ti = top_id[j]; top_id[j] = top_id[j+1]; top_id[j+1] = ti;
                ++j;
            }
        }
    }
    free(r_old);
    free(ids);
    free(ht_keys);
    free(ht_vals);

    // 7. 输出结果（从大到小）
    FILE* ofs = fopen("top100.txt", "w");
    if (!ofs) return 1;
    for (int i = K - 1; i >= 0; --i) {
        fprintf(ofs, "%u %.8f\n", top_id[i], top_val[i]);
        printf("%u %.8f\n", top_id[i], top_val[i]);
    }
    fclose(ofs);
    return 0;
}
