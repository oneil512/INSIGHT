import logging
import math
import os
import time
from collections import defaultdict, deque

import llama_index
from Bio import Entrez
from colorama import Fore
from langchain import OpenAI
from langchain.chat_models import ChatOpenAI
from llama_index import GPTListIndex, GPTSimpleVectorIndex, LLMPredictor, ServiceContext

from agents import boss_agent, worker_agent
from config import EMAIL, OPENAI_API_KEY
from utils import (
    execute_python,
    get_ada_embedding,
    get_code_params,
    get_key_results,
    insert_doc_llama_index,
    load,
    process_mygene_result,
    process_pubmed_result,
    query_knowledge_base,
    save,
    handle_python_result,
    handle_results
)

logging.getLogger("llama_index").setLevel(logging.WARNING)

Entrez.email = EMAIL or os.environ["EMAIL"]
MAX_TOKENS = 4097
api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]


def run(
    OBJECTIVE="",
    RESULT_CUTOFF=20000,
    MAX_ITERATIONS=1,
    TOOLS=["MYGENE", "PUBMED"],
    index=None,
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
                index,
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
        if not index:
            llm_predictor = LLMPredictor(
                llm=OpenAI(
                    temperature=0,
                    openai_api_key=api_key,
                    model_name="text-davinci-003",
                    max_tokens=2000,
                )
            )
            service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)
            index = GPTSimpleVectorIndex([], service_context=service_context)

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


    prev_iteration_no_result_notification = ""

    while True:
        # States used in each iteration
        result_code = None
        result_is_python = False

        task_list = boss_agent(
            objective=OBJECTIVE,
            tool_description=tool_description,
            task_list=task_list,
            index=index,
            completed_tasks=completed_tasks,
            no_result_notification=prev_iteration_no_result_notification
        )

        if not task_list:
            print(Fore.RED + "NO TASK LIST")
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
            result, cache, prev_iteration_no_result_notification = handle_python_result(result, cache)
            doc_store["tasks"][doc_store_task_key]["result_code"] = result_code

        if result:
            prev_iteration_no_result_notification = ""
            handle_results(result, index, doc_store, doc_store_task_key, task_id_counter, RESULT_CUTOFF)

        task_id_counter += 1

        if task_id_counter > MAX_ITERATIONS:
            break

    doc_store["key_results"] = get_key_results(index, OBJECTIVE, top_k=20)

    save(
        index,
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
run(OBJECTIVE=OBJECTIVE, MAX_ITERATIONS=MAX_ITERATIONS, TOOLS=TOOLS)

# Reload state and resume run
# TOOLS and MAX_ITERATIONS can also be passed in when reloading state.
# run(reload_path="out/Cure breast cancer_2023-04-25_16-38-42")
