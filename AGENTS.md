# FORBIDDEN reply/working style

- You are **FORBIDDEN** to tell what the user should do or verify: do and verify it yourself, and show it to user
- You are **FORBIDDEN** to answer **BEFORE** you actually run the code user provided or mentioned, and actually see what the error or output is.



# CRITICAL
## (CRITICAL) 用户交互规范
请用 **中文** 与用户交互
- **必须**: 你的**完整回答**必须**先**通过 **主对话窗口** 输出，**然后**调用 `question` 工具与用户交互
- **必须**: 只有当用户明确说明“任务完成”、“搞定”、“完事了”、“结束”等类似表达，才可结束
- **禁止**: 禁止主动结束对话或任务，只有当用户通过 `question` 工具明确回复"结束/完成/可以了"等意图时，才可结束
- **必须**: 当你即将完成一个请求或任务时，必须通过 `question` 工具请求用户反馈，确认是否满意或需要调整
- **必须**: 需求不明确时，必须通过 `question` 工具询问澄清，并提供预定义选项供用户选择
- **必须**: 存在多个方案/策略时，必须通过 `question` 工具列出选项让用户决定，禁止自作主张
- **必须**: 方案或策略需要变更时，必须通过 `question` 工具告知用户并获得确认
- **必须**: 必须预期用户有后续问题或者需求，请通过 `question` 工具与用户交互，放置预期选项，以及“任务完成，结束对话”选项。


### 调用 `question` 的函数
关键函数 `functions.question`

```bash
functions.question({
  "questions": [
    {
      "header": "下一步",
      "multiple": false,
      "options": [
        {
          "label": "继续找 null",
          "description": "继续计算缺失扇区的 null relations，逼出完整闭合"
        },
        {
          "label": "任务完成",
          "description": "目前结果够用，结束对话"
        }
      ],
      "question": "接下来你想让我做哪一步？"
    }
  ]
})
```

This is a non-optional protocol.

At every turn, you MUST perform exactly one of the following actions:

- Provide the requested output AND immediately call `question`, OR

- If any uncertainty exists, immediately call `question` without providing speculative output.

The conversation must never terminate voluntarily.

The assistant must never produce a terminal response.

`question` is mandatory at the end of every turn.




<!-- TRELLIS:START -->
# Trellis Instructions

These instructions are for AI assistants working in this project.

Use the `/trellis:start` command when starting a new session to:
- Initialize your developer identity
- Understand current project context
- Read relevant guidelines

Use `@/.trellis/` to learn:
- Development workflow (`workflow.md`)
- Project structure guidelines (`spec/`)
- Developer workspace (`workspace/`)

If you're using Codex, project-scoped helpers may also live in:
- `.agents/skills/` for reusable Trellis skills
- `.codex/agents/` for optional custom subagents

Keep this managed block so 'trellis update' can refresh the instructions.

<!-- TRELLIS:END -->
