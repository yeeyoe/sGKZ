# `k_stability_search` 数据库查看命令

下面的命令使用 SQLite 查看搜索数据库。先设置数据库路径：

```bash
DB=K-stability/k_stability_search.sqlite
```

如果使用自定义数据库路径，可以写成：

```bash
DB=./my-search.sqlite
```

## 数据库和表

```bash
ls -lh "$DB"
sqlite3 "$DB" '.tables'
sqlite3 "$DB" '.schema candidates'
sqlite3 "$DB" '.schema attempts'
sqlite3 "$DB" '.schema candidate_validations'
sqlite3 "$DB" '.schema state'

# 补算已有认证 witness 的精确 Q_P(g)^2（可重复执行）
./build/K-stability/k_stability_search --backfill-q --database "$DB"
```

## 候选数量

```bash
# candidates 总数
sqlite3 "$DB" 'select count(*) from candidates;'

# 指定 d 的候选总数
sqlite3 "$DB" \
  'select count(*) from candidates where d=3;'

# 按候选状态统计
sqlite3 -header -column "$DB" \
  'select status, count(*) as count
   from candidates
   group by status;'

# 指定 d 按状态统计
sqlite3 -header -column "$DB" \
  'select status, count(*) as count
   from candidates
   where d=3
   group by status;'

# verified_unstable 总数
sqlite3 "$DB" \
  'select count(*) from candidates
   where status="verified_unstable";'

# 指定 d 的 verified_unstable 总数
sqlite3 "$DB" \
  'select count(*) from candidates
   where d=3 and status="verified_unstable";'

# 指定 d 中被 Donaldson detector 测试过的不同候选总数
sqlite3 "$DB" \
  'select count(distinct a.candidate_key)
   from attempts as a join candidates as c on c.key=a.candidate_key
   where c.d=3;'

# 按 Q_P(g)^2 降序查看 verified 候选
sqlite3 -header -column "$DB" \
  'select key, twice_area, cast(twice_area as real)/2 as area
   from candidates
   where status="verified_unstable"
   order by q_squared_value desc, key;'

# 指定 d 按面积查看 verified 候选
sqlite3 -header -column "$DB" \
  'select key, twice_area, cast(twice_area as real)/2 as area
   from candidates
   where d=3 and status="verified_unstable"
   order by cast(twice_area as integer), key;'

# 指定 d 的 verified 候选前五名
sqlite3 -header -column "$DB" \
  'select key, twice_area
   from candidates
   where d=3 and status="verified_unstable"
   order by cast(twice_area as integer), key
   limit 5;'
```

## 指定候选

先设置候选的完整 key。例如：

```bash
KEY='d6|p=1:0;23:12;-19:10;-22:-5;-27:-29;91:-1|k=15,2,4,1,2,1'
```

查看候选的全部几何和精确数据：

```bash
sqlite3 -line "$DB" \
  "select key, d, directions, steps, vertices, normals,
          twice_area, status, ell0, ell1, ell2,
          q_squared_exact, q_squared_value
   from candidates
   where key='$KEY';"
```

字段含义：

- `directions`：方向序列 `p_1,...,p_d`，格式为 `x:y;...`；
- `steps`：步长序列 `k_1,...,k_d`，格式为 `k_1,k_2,...`；
- `vertices`：顶点序列，格式为 `x:y;...`；
- `normals`：facet normals 序列，格式为 `x:y;...`；
- `twice_area`：二倍面积，实际面积为 `twice_area/2`；
- `ell0,ell1,ell2`：`ell_P(x,y)=ell0+ell1*x+ell2*y`。
- `vertex_singularity_flags`：按顶点排列的 `0/1`，`0` 表示光滑，`1` 表示奇异；
- `singular_vertex_count`：奇异顶点总数，仅 verified 候选写入。
- `q_squared_exact`：候选所有认证 witness 中最大的精确 $Q_P(g)^2$；
- `q_squared_value`：同一指标的浮点表示，仅用于排序。

按奇异点数查看 verified 候选：

```bash
sqlite3 -header -column "$DB" \
  "select key, vertex_singularity_flags, singular_vertex_count
   from candidates
   where status='verified_unstable'
   order by singular_vertex_count, cast(twice_area as integer), key;"
```

## 检测和精确 witness

候选级状态决定搜索是否跳过，stage 级数值结果仍在 `attempts` 中：

```bash
sqlite3 -header -column "$DB" \
  "select candidate_key, validation_profile, status, last_stage, updated_at
   from candidate_validations
   where candidate_key='$KEY'
   order by updated_at;"
```

查看当前 profile 下尚未完成的候选：

```bash
sqlite3 -header -column "$DB" \
  "select candidate_key, status, last_stage
   from candidate_validations
   where validation_profile like 'validation-v2|%'
     and status='pending';"
```

查看该候选的所有 detector profile：

```bash
sqlite3 -header -column "$DB" \
  "select candidate_key, profile, status, value,
          witness_ux, witness_uy, witness_t,
          exact_a, exact_b, exact_c, exact_value, numerical_negative,
          q_squared_exact, q_squared_value
   from attempts
   where candidate_key='$KEY'
   order by rowid;"
```

只查看精确认证记录：

```bash
sqlite3 -header -column "$DB" \
  "select profile, exact_a, exact_b, exact_c, exact_value,
          q_squared_exact, q_squared_value
   from attempts
   where candidate_key='$KEY'
     and status='verified_unstable';"
```

只有 `verified_unstable` 记录应当含有 witness 字段；`unverified` 记录中的
witness 字段必须为 `NULL`。

attempts 中的 `q_squared_exact` 是该认证 witness 的精确有理数；候选表中的值
是所有认证 profile 中的最大值。

## 搜索断点状态

```bash
sqlite3 -header -column "$DB" \
  'select name, value from state order by name;'
```

其中包括按维数保存的随机数生成器状态，例如 `d3|rng`。shell 不作为断点保存，
每次运行都从 `shell=0` 开始；程序只读取当前 `d` 的候选和 detector profile
记录，不会加载或处理其它维数的候选。

## 结果报告

搜索 CLI 可以把当前最小 verified 候选写入指定目录：

```bash
./build/K-stability/k_stability_search \
  --d 6 --N 8 --M 8 --time-limit 3600 \
  --database "$DB" \
  --output-dir K-stability/K-results
```

报告文件为：

```text
K-stability/K-results/k_stability_search_result.txt
```

找到认证候选时，报告包含完整几何、面积、边界测度、`ell_P`、检测 profile
和精确 witness；没有找到时文件内容为：

```text
没找到
```
