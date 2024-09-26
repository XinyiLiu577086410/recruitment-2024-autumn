# 七边形 2024 年春季招新题目！

嗨，欢迎来做七边形的招新题！

招新题是一个C++ 小项目，使用 CMake 构建系统。你需要对这个项目做你能想到的性能优化，使其在保持正确性的前提下，跑得尽可能快。

* 我们将会提供七边形的机器供大家在上面编程和调试。最终加速比（也就是性能表现）将会以在七边形机器上测试的结果为准。
* 集群使用 [spack](https://spack.io/) 进行包管理，你可以通过spack找到一些你想要的包。
* 默认的机器（即登录节点）仅有一台，且没有 GPU。
* 但队内集群中数个节点包含GPU。集群通过 [slurm](https://slurm.schedmd.com/documentation.html) 进行管理，如果你编写了很酷的多机并行算法，或者编写了很酷的 GPU 加速模块，请使用 slurm 运行你的应用程序。
* 如果遇到无法解决的问题，可以联系管理员获取相关支持。

请仔细阅读下面的规则、题目详情和提交方式。

截止时间为 **2024年10月14号**。你需要要那之前完成并提交所有工作。


## 题目大意
**序列比对(Sequence Alignment)**是生物信息学中最重要的技术之一。序列比对的结果是许多其他步骤的基础。它可以用来找出对齐序列之间的差异和相似性，这是生物序列识别、结构预测和功能分析的前提。然而，序列比对是一项非常耗时的任务。近年来，随着并行比对算法的不断成熟和优化，计算时间显著减少。

根据所使用的比对方法，序列比对算法通常可以分为两种主要类型：全局比对和局部比对。1970年，Saul B. Needleman和Christian D. Wunsch提出了**Needleman-Wunsch(NW)算法**，以找到整个序列中的最佳匹配。随后，在1981年，Temple F. Smith和Michael S. Waterman基于NW算法开发了局部比对算法，后来称为**Smith-Waterman(SW)算法**，用于寻找序列对之间的最佳子序列匹配。NW和SW算法都应用动态规划(DP)来计算序列比对，这使得这两种算法的时间复杂度为平方级。因此，大规模序列比对在计算上非常消耗资源。为了加速这一过程，学术界/工业界进行了大量努力。大多数任务通过并行计算加速，包括向量级、线程级、进程级和异构并行化。

本次招新的笔试任务是对**Smith-Waterman算法**进行优化。该算法实现了局部的序列比对（local sequence alignment），即确定两个核酸序列或蛋白质序列之间的相似区域。具体而言，本次的算例为将人类第17号染色体中的部分核酸序列与五个物种（倭黑猩猩（Bonobo，学名 Pan paniscus）、黑猩猩（Chimpanzee，学名 Pan troglodytes）、大猩猩（Gorilla）、东非狒狒（Olive baboon，学名 Papio anubis）和苏门答腊猩猩（Sumatran orangutan，学名 Pongo abelii））的第17号染色体部分核酸序列进行局部比对。由于完整染色体的核酸序列过于庞大，我们仅截取了这些物种的不同规模的一小部分序列，因此最终的比对结果不一定与实际情况相符。

## 实现要求

项目目录如下：

```
.
├── CMakeLists.txt
├── data
├── include
│   ├── SmithWaterman.hpp
│   └── Timer.hpp
├── main.cc
├── README.md
└── src
    └── SmithWaterman.cc
```

* `data`目录存放了计算所需的输入数据及参考文件（禁止更改该目录下的所有文件）。我们目前提供了3个规模的算例，baseline的时间分别为`8.72s`, `35.32s`, `3880.45s`。之后可能会增加新的算例。
* 在**不更改计算逻辑**（即原来有的步骤不能直接删掉）且**结果正确**（data目录下的参考文件禁止更改）的情况下，可以任意修改项目源文件及构建文件。


### 注意
* 本题目主要考察对该算法的优化，故**禁止**调用任何已经实现Smith-Waterman算法的库。
* **禁止**根据输入数据的某些特征对程序进行优化，例如打表等。
* **禁止**更改`match_score`, `mismatch_score` 和 `gap_score` 的值。




## 编译、运行和测试

通过
```
cmake -B build && cmake --build build
```
构建项目。

通过
```
ctest --test-dir build
```
直接在命令行执行。

但有些算例的执行时间较长。~~若是在执行期间突然断连，就会发生可怕的事情。~~因此，无论是否通过多节点/GPU执行你的程序，我们都建议通过slurm提交你的任务（也方便我们对你的程序进行测试）。


最终评测时，我们会通过
```
cmake -B build && cmake --build build
```
构建你的程序。

然后通过
```
sbatch job.slurm
```
将你的任务提交给slurm执行。



## 提交

你需要 fork 题目仓库，并在代码提交截止日期前将自己的仓库链接以 issue 的方式发表在本仓库下。你的仓库应该包含优化实现的代码。

与此同时，你需要制作一份 ppt，内容如下：

* 每一步优化的思路、过程和效果（举例：使用了 xxx 优化，相对于原代码，速度提升了 114.514 倍）
* 最好对程序进行 profile，以了解性能瓶颈
* 你在解题过程中所参考的资料（如有使用人工智能工具，请注明）
* 在解题过程中，遇到的有意思的事情，或者是让你印象深刻的 bug（可选）

代码审查会用到这份 ppt，因此请你把它发送到超算队邮箱。

另外，在笔试结束后（具体时间另外通知），你需要准备一次展示，向我们介绍你的优化成果。

如果对题目本身或者提交方式有任何问题，请积极在群里讨论。

## 提示

* 良好的版本控制习惯可以有效避免 “改了一下午代码后它跑不起来了，但是我又改不回去了” 之类的极端情况。推荐使用 [Git](https://git-scm.com/)
* 使用 htop/glances 或者资源管理器查看 CPU 负载是验证自己程序正在并行工作的好方法之一
* Perf 是 Linux 下常用的性能分析工具：<https://www.brendangregg.com/perf.html>
* 火焰图可以直观地展示程序的性能热点：<https://www.brendangregg.com/flamegraphs.html>
* OpenMP 是常用的并行计算框架：<https://bisqwit.iki.fi/story/howto/openmp/>
* 如果需要做多机的并行，你可能会对 MPI 感兴趣：<https://mpitutorial.com/tutorials/>
* SIMD 是提升程序单核性能的有效手段：<https://zhuanlan.zhihu.com/p/583326378> <https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#>
* 适量的 Cache 知识不仅对单核优化有用，也能避免一些多线程场景下的诡异问题：<https://zhuanlan.zhihu.com/p/136300660>
* 更加系统性的高性能计算知识，请参考七边形的 HPC 学习路线：<https://heptagonhust.github.io/HPC-roadmap/>
* Slurm 怎么用捏 <http://scc.ustc.edu.cn/zlsc/user_doc/html/slurm/index.html>
* 什么是 Smith Waterman 算法捏 <https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm>
* 怎么并行 Smith Waterman 算法捏 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8419822/#CR37>


以及一些技术无关的提示：

* 如果遇到问题，请积极在 QQ 群内提问
* 部分优化可能相当难做/难调。所以即使没有写出程序/调试成功，也欢迎把自己的天才优化想法写进提交的文档里，通过文字来展示自己对高性能计算的理解和认知。
