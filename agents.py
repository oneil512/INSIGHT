import os
from ast import literal_eval
from collections import deque
from typing import List

import openai

from utils import generate_tool_prompt, get_gpt_chat_completion, get_gpt_completion

openai.api_key = os.environ["OPENAI_API_KEY"]
openai.organization = os.environ["OPENAI_ORG"]
tools = ["MYGENE", "PUBMED"]


def boss_agent(
    objective: str,
    tool_description: str,
    task_list: List[str],
    executive_summary="No tasks completed Yet.",
    completed_tasks=List[str],
):
    system_prompt = """You are BossGPT, a responsible and organized agent that is responsible for completing a high level and difficult objective. 
As the boss, your goal is to break the high level objective down into small and managable tasks for your workers. These tasks will be picked up by your worker agents and completed. 
You will also get an executive summary of what your workers have accomplished so far. Use the summary to make decisions about what tasks to do next, what tasks to get rid of, 
and to reprioritize tasks. The highest priority task will be at the top of the task list.

You also have access to some tools. You can create a task for your workers to use any of your tools. You cannot use more than one tool per task.
Your worker agents update the executive summary so that you can use new information from the completed tasks to make informed decisions about what to do next. 
It is ok to create tasks that do not directly help achieve your objective but rather just serve to add useful information.

Tasks should be a simple python array with strings as elements. Priority is only determined by the order of the elements in the array.

After you have finished generating your task array, cease output.

===

Your responses should in this format:

THOUGHTS
{Reason about what tasks to add, change, delete, or reprioritize given your objective and all the information you have}

TASKS
{Python array of tasks}

""".strip()

    user_prompt = f"""
Here is your objective: {objective}
Here is the current task list: {task_list}
Here are the tasks that have been complete thus far: {completed_tasks}
Here are the tools you have access to: {tool_description}
Here is an executive summary of the information gathered so far {executive_summary}

If a task has already been completed, do not write that same task again in the task list. If you would like a worker to continue or redo a task, be sure to word it a little differently so you don't get the same result.

===

Please update the task list and follow this format.

THOUGHTS
Reason about what tasks to add, change, delete, or reprioritize given your objective and the information you have

TASKS
Python array of tasks

===

Here is an example of the tasks list. Be sure that it is valid python:

TASKS
["Research frog habitats", "Find all species of trees", "Get world population", "Retrieve facts about the american civil war"]

Note: To be sure that TASKS is a valid python list, it should always start with '[' and always end with ']'
    
    """.strip()

    content = get_gpt_chat_completion(system_prompt, user_prompt, temp=0.0)

    thoughts = content[
        content.find("THOUGHTS") + len("THOUGHTS") : content.find("TASKS")
    ].strip()
    tasks = content[content.find("TASKS") + len("TASKS") :].strip()

    # parsed_tasks = parser("Parse the following text so that it is a valid python list. Do not alter the elements in any way.", tasks)
    new_task_list = literal_eval(tasks)
    return deque(new_task_list), thoughts


def worker_agent(
    objective: str,
    task: str,
    context: str,
    previous_params: str = None,
    python: bool = False,
) -> str:
    if python:
        prompt = f"""You are an AI who performs one task based on the following objective: {objective}. You will be writing code for your task. Here are the parameters used in the code from previous tasks {previous_params}. Do not use the same parameters and query again; instead use tweak the parameters or query so that you get a slightly different result."""
        prompt += generate_tool_prompt(task)
    else:
        prompt = f"""You are an AI who performs one task based on the following objective: {objective}. Here is the result of some similar previous tasks: {context}. Try to not produce the same result. Be creative so that we can get new information."""

    prompt += f"\nYour task: {task}\nResponse:"

    response = get_gpt_completion(prompt, engine="text-davinci-003", temp=0.0)

    return response


def data_cleaning_agent(result: str, objective: str) -> str:
    prompt = f"""You are an AI who summarizes, cleans, and organizes data. It is important that you do not delete any information that could be useful. Respond with only the updated information.
Data: {result}
Cleaned Data:"""
    response = get_gpt_completion(prompt, engine="text-davinci-003", temp=0.1)

    return response
