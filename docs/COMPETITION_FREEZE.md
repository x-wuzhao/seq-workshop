# 参赛冻结与独立迭代边界

用户于2026-11-01明确要求：比赛已经进入省赛且提交内容锁定；后续开发作为独立产品保存在 GitHub，不改变比赛链接及其内容。

- 原参赛仓库：x-wuzhao/seq-workshop
- 保留分支：main
- 保留提交：57e9a4f2ba14234026f9417019038509c6f6c3de
- 比赛地址：<https://x-wuzhao.github.io/seq-workshop/>
- 迭代分支：product-next，从保留提交分出，纳入本地后续开发并继续修正。
- 本轮开始时页面 UTF-8 内容 SHA-256：ceb75ee0bacfe528862b5e97c5c8ec2e7e745d69298bde07188d93710029190d。

本分支不合并 main，不修改 Pages 设置，不使用 gh-pages 发布，不配置 Pages deploy workflow。CI 只在 product-next 的 push/PR 运行，权限仅 contents:read，输出离线构建 artifact。

之后如需发布新产品，应使用独立站点或仓库，并先明确新地址；不得复用比赛网页的部署入口。即使 URL 字符串不变，更新网页内容也属于改变参赛版本。

原始本地源码目录 D:\ChatGpt\seq-workshop 保留；迭代检出为 D:\ChatGpt\seq-workshop-next。工作区备份键与 Service Worker 缓存前缀和旧版分开。
