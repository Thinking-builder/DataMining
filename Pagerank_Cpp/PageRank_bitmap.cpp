// // #include <cstdio>
// // #include <cstdlib>
// // #include <cstdint>
// // #include <cmath>
// // #include <cstring>

// // const double BETA = 0.85;
// // const double EPSILON = 1e-4;

// // int main() {
// //     // 1. 找出最大ID
// //     FILE *fin = fopen("Data.txt", "r");
// //     if (!fin) return 1;
// //     int u, v;
// //     int max_id = 0;
// //     while (fscanf(fin, "%d %d", &u, &v) == 2) {
// //         if (u > max_id) max_id = u;
// //         if (v > max_id) max_id = v;
// //     }
// //     // 位图大小
// //     int bits = (max_id >> 3) + 1;
// //     uint8_t *bitset = (uint8_t*)calloc(bits, 1);
// //     // 2. 收集唯一ID
// //     rewind(fin);
// //     int N = 0;
// //     int *id_list = (int*)malloc((max_id + 1) * sizeof(int)); // 临时大小
// //     while (fscanf(fin, "%d %d", &u, &v) == 2) {
// //         if (u >= 0) {
// //             int idx = u >> 3, off = u & 7;
// //             if (!(bitset[idx] & (1 << off))) {
// //                 bitset[idx] |= (1 << off);
// //                 id_list[N++] = u;
// //             }
// //         }
// //         if (v >= 0) {
// //             int idx = v >> 3, off = v & 7;
// //             if (!(bitset[idx] & (1 << off))) {
// //                 bitset[idx] |= (1 << off);
// //                 id_list[N++] = v;
// //             }
// //         }
// //     }
// //     // 排序并唯一（可能重复添加）
// //     fclose(fin);
// //     // 简单排序
// //     for (int i = 0; i < N - 1; ++i) {
// //         for (int j = i + 1; j < N; ++j) {
// //             if (id_list[j] < id_list[i]) {
// //                 int t = id_list[i]; id_list[i] = id_list[j]; id_list[j] = t;
// //             }
// //         }
// //     }
// //     int uniq = 0;
// //     for (int i = 0; i < N; ++i) {
// //         if (i == 0 || id_list[i] != id_list[i - 1]) {
// //             id_list[uniq++] = id_list[i];
// //         }
// //     }
// //     N = uniq;
// //     free(bitset);

// //     // 3. 保存索引至id的映射
// //     int *idx2id = id_list; // 0..N-1 maps to original id

// //     // 分配PageRank和degree
// //     uint16_t *degree = (uint16_t*)calloc(N, sizeof(uint16_t));
// //     float *r_old = (float*)malloc(N * sizeof(float));
// //     float *r_new = (float*)malloc(N * sizeof(float));
// //     if (!degree || !r_old || !r_new) return 1;

// //     // 二分查找函数
// //     auto map_id = [&](int x) {
// //         int l = 0, r = N - 1;
// //         while (l <= r) {
// //             int m = (l + r) >> 1;
// //             if (idx2id[m] == x) return m;
// //             if (idx2id[m] < x) l = m + 1;
// //             else r = m - 1;
// //         }
// //         return -1;
// //     };

// //     // 4. 第三遍扫描：统计出度
// //     fin = fopen("Data.txt", "r");
// //     while (fscanf(fin, "%d %d", &u, &v) == 2) {
// //         int ui = map_id(u);
// //         if (ui >= 0) degree[ui]++;
// //     }

// //     // 初始PageRank
// //     for (int i = 0; i < N; ++i) r_old[i] = 1.0f / N;

// //     // PageRank迭代（重新扫描文件）
// //     while (1) {
// //         for (int i = 0; i < N; ++i) r_new[i] = (1.0f - BETA) / N;
// //         float dsum = 0;
// //         for (int i = 0; i < N; ++i) if (degree[i] == 0) dsum += r_old[i];
// //         rewind(fin);
// //         while (fscanf(fin, "%d %d", &u, &v) == 2) {
// //             int ui = map_id(u), vi = map_id(v);
// //             if (ui >= 0 && vi >= 0 && degree[ui]) {
// //                 r_new[vi] += BETA * (r_old[ui] / degree[ui]);
// //             }
// //         }
// //         float avg = BETA * dsum / N;
// //         for (int i = 0; i < N; ++i) r_new[i] += avg;
// //         float diff = 0;
// //         for (int i = 0; i < N; ++i) diff += fabsf(r_new[i] - r_old[i]);
// //         if (diff < EPSILON) break;
// //         memcpy(r_old, r_new, N * sizeof(float));
// //     }
// //     fclose(fin);

// //     // Top-100 插入排序
// //     int top_n = N < 100 ? N : 100;
// //     int *top_id = (int*)malloc(top_n * sizeof(int));
// //     float *top_val = (float*)malloc(top_n * sizeof(float));
// //     int count = 0;
// //     for (int i = 0; i < N; ++i) {
// //         float v = r_old[i];
// //         int k = count;
// //         if (count < top_n) {
// //             top_id[count] = i;
// //             top_val[count++] = v;
// //         } else if (v > top_val[0]) {
// //             top_val[0] = v;
// //             top_id[0] = i;
// //             k = 0;
// //         } else continue;
// //         while (k + 1 < count && top_val[k] > top_val[k + 1]) {
// //             float tv = top_val[k+1]; top_val[k+1] = top_val[k]; top_val[k] = tv;
// //             int ti = top_id[k+1]; top_id[k+1] = top_id[k]; top_id[k] = ti;
// //             ++k;
// //         }
// //     }
// //     // 输出前100
// //     FILE *ofs = fopen("top100.txt","w");
// //     for (int i = count-1; i >= 0; --i) {
// //         int oid = idx2id[top_id[i]];
// //         fprintf(ofs, "Node %d: %.8f\n", oid, top_val[i]);
// //         printf("Node %d: %.8f\n", oid, top_val[i]);
// //     }
// //     fclose(ofs);
// //     return 0;
// // }

// #include <cstdio>
// #include <cstdlib>
// #include <cstdint>
// #include <cmath>

// // qsort 比较函数
// static int cmp_uint32(const void* a, const void* b) {
//     uint32_t A = *(const uint32_t*)a;
//     uint32_t B = *(const uint32_t*)b;
//     return (A > B) - (A < B);
// }

// // 二分查找
// static int bsearch_uint32(const uint32_t* arr, int n, uint32_t key) {
//     int l = 0, r = n - 1;
//     while (l <= r) {
//         int m = (l + r) >> 1;
//         if (arr[m] == key) return m;
//         if (arr[m] < key) l = m + 1;
//         else r = m - 1;
//     }
//     return -1;
// }

// static const double BETA = 0.85;
// static const double EPS  = 1e-4;

// int main() {
//     FILE* fin = fopen("Data.txt", "r");
//     if (!fin) return 1;

//     // 1. 读所有端点到动态数组
//     int cap = 1 << 18; // 初始容量 262144
//     int cnt = 0;
//     uint32_t *ids = (uint32_t*)malloc(sizeof(uint32_t) * cap);
//     uint32_t u, v;
//     while (fscanf(fin, "%u %u", &u, &v) == 2) {
//         if (cnt + 2 > cap) {
//             cap <<= 1;
//             ids = (uint32_t*)realloc(ids, sizeof(uint32_t) * cap);
//         }
//         ids[cnt++] = u;
//         ids[cnt++] = v;
//     }

//     // 2. 排序去重
//     qsort(ids, cnt, sizeof(uint32_t), cmp_uint32);
//     int N = 0;
//     for (int i = 0; i < cnt; ++i) {
//         if (i == 0 || ids[i] != ids[i-1]) {
//             ids[N++] = ids[i];
//         }
//     }
//     ids = (uint32_t*)realloc(ids, sizeof(uint32_t) * N);

//     // 3. 统计出度 (uint16_t)
//     uint16_t *degree = (uint16_t*)calloc(N, sizeof(uint16_t));
//     rewind(fin);
//     while (fscanf(fin, "%u %u", &u, &v) == 2) {
//         int ui = bsearch_uint32(ids, N, u);
//         if (ui >= 0) degree[ui]++;
//     }

//     // 4. 分配 PR 数组
//     float *r_old = (float*)malloc(sizeof(float) * N);
//     float *r_new = (float*)malloc(sizeof(float) * N);
//     for (int i = 0; i < N; ++i) r_old[i] = 1.0f / N;

//     // 5. 迭代计算 PageRank
//     while (1) {
//         double dsum = 0;
//         for (int i = 0; i < N; ++i) {
//             if (degree[i] == 0) dsum += r_old[i];
//             r_new[i] = (1.0 - BETA) / N;
//         }
//         rewind(fin);
//         while (fscanf(fin, "%u %u", &u, &v) == 2) {
//             int ui = bsearch_uint32(ids, N, u);
//             if (ui < 0 || degree[ui] == 0) continue;
//             int vi = bsearch_uint32(ids, N, v);
//             r_new[vi] += (float)(BETA * (r_old[ui] / degree[ui]));
//         }
//         float avg = (float)(BETA * dsum / N);
//         double diff = 0;
//         for (int i = 0; i < N; ++i) {
//             r_new[i] += avg;
//             diff += fabs(r_new[i] - r_old[i]);
//         }
//         if (diff < EPS) break;
//         float *tmp = r_old;
//         r_old = r_new;
//         r_new = tmp;
//     }
//     fclose(fin);
//     free(degree);
//     free(r_new);

//     // 6. 在线保持 Top-100
//     int K = N < 100 ? N : 100;
//     float  top_val[100];
//     uint32_t top_id[100];
//     int cnt_top = 0;
//     for (int i = 0; i < N; ++i) {
//         float val = r_old[i];
//         uint32_t id = ids[i];
//         if (cnt_top < K) {
//             int j = cnt_top++;
//             top_val[j] = val;
//             top_id[j]  = id;
//             while (j > 0 && top_val[j-1] > top_val[j]) {
//                 // 小值上移
//                 float tv = top_val[j]; top_val[j] = top_val[j-1]; top_val[j-1] = tv;
//                 uint32_t ti = top_id[j]; top_id[j] = top_id[j-1]; top_id[j-1] = ti;
//                 --j;
//             }
//         } else if (val > top_val[0]) {
//             top_val[0] = val;
//             top_id[0]  = id;
//             int j = 0;
//             while (j + 1 < K && top_val[j] > top_val[j+1]) {
//                 float tv = top_val[j]; top_val[j] = top_val[j+1]; top_val[j+1] = tv;
//                 uint32_t ti = top_id[j]; top_id[j] = top_id[j+1]; top_id[j+1] = ti;
//                 ++j;
//             }
//         }
//     }
//     free(r_old);
//     free(ids);

//     // 7. 输出 Top-100（从大到小）
//     FILE* ofs = fopen("top100.txt", "w");
//     if (!ofs) return 1;
//     for (int i = K - 1; i >= 0; --i) {
//         fprintf(ofs, "Node %u: %.8f\n", top_id[i], top_val[i]);
//         printf("Node %u: %.8f\n", top_id[i], top_val[i]);
//     }
//     fclose(ofs);
//     return 0;
// }


#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>

static const double BETA = 0.85;
static const double EPS  = 1e-4;

// 位图宏
#define BIT_GET(bs, i)   ((bs[(i) >> 3] >> ((i) & 7)) & 1)
#define BIT_SET(bs, i)   (bs[(i) >> 3] |= (1 << ((i) & 7)))

// 快速比较函数，用于 qsort
static int cmp_u32(const void* a, const void* b) {
    uint32_t A = *(const uint32_t*)a;
    uint32_t B = *(const uint32_t*)b;
    return (A > B) - (A < B);
}

// 二分查找返回索引或 -1
static int find_idx(const uint32_t* arr, int n, uint32_t key) {
    int l = 0, r = n - 1;
    while (l <= r) {
        int m = (l + r) >> 1;
        if (arr[m] == key) return m;
        if (arr[m] < key) l = m + 1;
        else r = m - 1;
    }
    return -1;
}

int main() {
    // 打开文件
    FILE* fin = fopen("Data.txt", "r");
    if (!fin) return 1;
    uint32_t u, v;
    uint32_t max_id = 0;
    size_t M = 0;
    // 第1遍：找 max_id 和边数
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        if (u > max_id) max_id = u;
        if (v > max_id) max_id = v;
        ++M;
    }
    // 位图大小
    uint32_t bits = max_id + 1;
    uint32_t bytes = (bits + 7) >> 3;
    uint8_t* bitset = (uint8_t*)calloc(bytes, 1);
    if (!bitset) return 1;
    // 第2遍：标记唯一节点
    rewind(fin);
    int N = 0;
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        if (!BIT_GET(bitset, u)) { BIT_SET(bitset, u); ++N; }
        if (!BIT_GET(bitset, v)) { BIT_SET(bitset, v); ++N; }
    }
    // 构建 ids
    uint32_t* ids = (uint32_t*)malloc(N * sizeof(uint32_t));
    int idx = 0;
    for (uint32_t i = 0; i < bits; ++i) {
        if (BIT_GET(bitset, i)) ids[idx++] = i;
    }
    free(bitset);
    // 排序 ids
    qsort(ids, N, sizeof(uint32_t), cmp_u32);

    // 第3遍：统计出度
    uint32_t* degree = (uint32_t*)calloc(N, sizeof(uint32_t));
    if (!degree) return 1;
    rewind(fin);
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        int ui = find_idx(ids, N, u);
        if (ui >= 0) degree[ui]++;
    }
    fclose(fin);

    // 构建 CSR 偏移
    uint32_t* offset = (uint32_t*)malloc((N + 1) * sizeof(uint32_t));
    offset[0] = 0;
    for (int i = 1; i <= N; ++i) {
        offset[i] = offset[i-1] + degree[i-1];
    }
    uint32_t E = offset[N];
    // 邻接表
    uint16_t* adj = (uint16_t*)malloc(E * sizeof(uint16_t));
    // 重用 degree 作为 head
    memcpy(degree, offset, N * sizeof(uint32_t));
    // 第4遍：填充邻接表
    fin = fopen("Data.txt", "r");
    if (!fin) return 1;
    while (fscanf(fin, "%u %u", &u, &v) == 2) {
        int ui = find_idx(ids, N, u);
        int vi = find_idx(ids, N, v);
        if (ui >= 0 && vi >= 0) {
            uint32_t pos = degree[ui]++;
            adj[pos] = (uint16_t)vi;
        }
    }
    fclose(fin);
    free(degree);

    // 初始化 PageRank
    float* r_old = (float*)malloc(N * sizeof(float));
    float* r_new = (float*)malloc(N * sizeof(float));
    float invN = 1.0f / N;
    for (int i = 0; i < N; ++i) r_old[i] = invN;

    // 迭代 PageRank
    while (1) {
        double dsum = 0.0;
        for (int i = 0; i < N; ++i) {
            uint32_t deg_i = offset[i+1] - offset[i];
            if (deg_i == 0) dsum += r_old[i];
            r_new[i] = (1.0f - BETA) * invN;
        }
        for (int i = 0; i < N; ++i) {
            uint32_t deg_i = offset[i+1] - offset[i];
            if (deg_i == 0) continue;
            float contrib = (float)(BETA * (r_old[i] / deg_i));
            for (uint32_t k = offset[i]; k < offset[i+1]; ++k) {
                r_new[adj[k]] += contrib;
            }
        }
        float avg = (float)(BETA * dsum * invN);
        double diff = 0.0;
        for (int i = 0; i < N; ++i) {
            r_new[i] += avg;
            diff += fabsf(r_new[i] - r_old[i]);
        }
        if (diff < EPS) break;
        float* tmp = r_old; r_old = r_new; r_new = tmp;
    }
    free(r_new);
    free(offset);

    // Top-100 在线选择
    int K = N < 100 ? N : 100;
    float top_val[100];
    uint32_t top_id[100];
    int cnt_top = 0;
    for (int i = 0; i < N; ++i) {
        float val = r_old[i];
        uint32_t idv = ids[i];
        if (cnt_top < K) {
            int j = cnt_top++;
            top_val[j] = val;
            top_id[j] = idv;
            while (j > 0 && top_val[j-1] > top_val[j]) {
                float tv = top_val[j]; top_val[j] = top_val[j-1]; top_val[j-1] = tv;
                uint32_t ti = top_id[j]; top_id[j] = top_id[j-1]; top_id[j-1] = ti;
                --j;
            }
        } else if (val > top_val[0]) {
            top_val[0] = val;
            top_id[0]  = idv;
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
    free(adj);

    // 输出结果
    FILE* ofs = fopen("top100.txt", "w");
    if (!ofs) return 1;
    for (int i = K - 1; i >= 0; --i) {
        fprintf(ofs, "Node: %u %.8f\n", top_id[i], top_val[i]);
        printf("%u %.8f\n", top_id[i], top_val[i]);
    }
    fclose(ofs);
    return 0;
}
