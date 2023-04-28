import logging
import math
import os
import time
from collections import defaultdict, deque

from Bio import Entrez
from colorama import Fore

from agents import boss_agent, worker_agent
from config import EMAIL, OPENAI_API_KEY
from utils import (
    get_key_results,
    insert_doc_llama_index,
    load,
    query_knowledge_base,
    save,
    handle_python_result,
    handle_results,
    create_index,
)

logging.getLogger("llama_index").setLevel(logging.WARNING)

Entrez.email = EMAIL or os.environ["EMAIL"]
MAX_TOKENS = 4097
api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]

indicies = []
summaries = []


def run(
    OBJECTIVE="",
    RESULT_CUTOFF=20000,
    MAX_ITERATIONS=1,
    TOOLS=["MYGENE", "PUBMED"],
    master_index=None,
    task_id_counter=1,
    task_list=deque(),
    completed_tasks=[],
    cache=defaultdict(list),
    current_datetime="",
    reload_path="",
):
    start_time = time.time()
    reload_count = 0
    if reload_path:
        try:
            (
                master_index,
                task_id_counter,
                task_list,
                completed_tasks,
                cache,
                current_datetime,
                OBJECTIVE,
                reload_count,
            ) = load(reload_path)
            MAX_ITERATIONS += task_id_counter - 1
        except Exception as e:
            print(f"Issue loading state {e}, with path {reload_path}")
            raise Exception("Cannot reload state.")

    else:
        if not master_index:
            master_index = create_index(api_key=api_key)
        current_datetime = current_datetime or str(time.strftime("%Y-%m-%d_%H-%M-%S"))

    doc_store = {"tasks": {}}

    tool_description_mapping = {
        "PUBMED": """2) Query PubMed API. This is useful for searching biomedical literature and studies on any medical subject. If you wish to make a task to create an API request to the PubMed API then simply say 'PUBMED:' followed by what you would like to search for. Example: 'PUBMED: Find recent developments in HIV research'""",
        "MYGENE": """1) Query mygene API. This is useful for finding information on specific genes, or genes associated with the search query. If you wish to make a task to create an API request to mygene then simply say 'MYGENE:' followed by what you would like to search for. Example: 'MYGENE: look up information on genes that are linked to cancer'""",
    }

    tool_description = ""

    for tool in TOOLS:
        tool_description += f"{tool_description_mapping[tool]}\n"

    print(Fore.CYAN, "\n*****OBJECTIVE*****\n")
    print(OBJECTIVE)


    result = []
    task = ""
    index = create_index(api_key)

    while True:
        result_is_python = False

        task_list = boss_agent(
            objective=OBJECTIVE,
            tool_description=tool_description,
            task_list=task_list,
            summaries=summaries,
            completed_tasks=completed_tasks,
            previous_task=task,
            previous_result=result
        )


        if not task_list:
            print(Fore.RED + "TASK LIST EMPTY")
            break

        print(Fore.WHITE + "\n*****TASK LIST*****\n")
        for t in task_list:
            print(t)

        task = task_list.popleft()

        print(Fore.RED + "\n*****NEXT TASK*****\n")
        print("task id: ", task_id_counter, "task: ", task)

        if any(tool in task for tool in TOOLS):
            result_is_python = True

        result = worker_agent(OBJECTIVE, task, index, cache, TOOLS, result_is_python)
        completed_tasks.append(task)

        print(Fore.GREEN + "\n*****TASK RESULT*****\n")
        print(Fore.GREEN + result)

        doc_store_task_key = str(task_id_counter) + "_" + task
        doc_store["tasks"][doc_store_task_key] = {}
        doc_store["tasks"][doc_store_task_key]["results"] = []

        if result_is_python:
            result = handle_python_result(result, cache, task, doc_store, doc_store_task_key)

        handle_results(result, index, doc_store, doc_store_task_key, task_id_counter, RESULT_CUTOFF)

        if index.docstore.docs:
            query="Give a detailed but terse overview of all the information. Start with a high level summary and then go into details. Do not include any further instruction. Do not include filler words. Include citation information."
            executive_summary = query_knowledge_base(index, query=query, list_index=False)
            insert_doc_llama_index(index=master_index, doc_id=str(task_id_counter), metadata=executive_summary)
            summaries.append(executive_summary)
            indicies.append(index)
            doc_store["tasks"][doc_store_task_key]["executive_summary"] = executive_summary

            index = create_index(api_key=api_key)

        task_id_counter += 1

        if task_id_counter > MAX_ITERATIONS:
            break

    doc_store["key_results"] = get_key_results(master_index, OBJECTIVE, top_k=1)

    save(
        master_index,
        doc_store,
        OBJECTIVE,
        current_datetime,
        task_id_counter,
        task_list,
        completed_tasks,
        cache,
        reload_count,
    )

    end_time = time.time()
    total_time = end_time - start_time
    print(Fore.RED + f"Total run time: {total_time:.2f} seconds")



### Set variables here.

TOOLS = ["MYGENE", "PUBMED"]
MAX_ITERATIONS = 1
OBJECTIVE = "Cure breast cancer"


# If you would like to reload a previous state, comment out run(OBJECTIVE=OBJECTIVE, MAX_ITERATIONS=MAX_ITERATIONS, TOOLS=TOOLS) and uncomment #run(reload_path="out/Cure breast cancer_2023-04-25_16-38-42")
# Then put your path in to your saved state.
# Note that state reloading is not backwards compatible. Saved states before this change cannot be reloaded.

# Fresh Run
#run(OBJECTIVE=OBJECTIVE, MAX_ITERATIONS=MAX_ITERATIONS, TOOLS=TOOLS)

# Reload state and resume run
# TOOLS and MAX_ITERATIONS can also be passed in when reloading state.
run(reload_path="out/Cure breast cancer_2023-04-28_01-01-23")
