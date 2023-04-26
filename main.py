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
        python = False
        cache_params = None
        params = None

        if index.docstore.docs:
            executive_summary = query_knowledge_base(index)
        else:
            executive_summary = "No information gathered yet."

        task_list, thoughts = boss_agent(
            objective=OBJECTIVE,
            tool_description=tool_description,
            task_list=task_list,
            executive_summary=executive_summary,
            completed_tasks=completed_tasks,
            no_result_notification=prev_iteration_no_result_notification
        )

        print(Fore.CYAN + "\n*****EXECUTIVE SUMMARY*****\n")
        print(Fore.CYAN + executive_summary)

        print(Fore.CYAN + "\n*****BOSS THOUGHTS*****\n")
        print(Fore.CYAN + thoughts)

        if task_list:
            print(Fore.WHITE + "\n*****TASK LIST*****\n")

            for t in task_list:
                print(t)

            task = task_list.popleft()

            context = ""
            if index.docstore.docs:
                context = query_knowledge_base(
                    index,
                    query=f"Provide as much useful context as possible for this task: {task}",
                )

            if any(tool in task for tool in TOOLS):
                python = True

            if "MYGENE" in task and "MYGENE" in cache:
                cache_params = cache["MYGENE"]

            if "PUBMED" in task and "PUBMED" in cache:
                cache_params = cache["PUBMED"]

            print(Fore.RED + "\n*****NEXT TASK*****\n")
            print("task id: ", task_id_counter, "task: ", task)

            result = worker_agent(OBJECTIVE, task, context, cache_params, python)
            completed_tasks.append(task)
            print(Fore.GREEN + "\n*****TASK RESULT*****\n")
            print(Fore.GREEN + result)

            # If task result was python code
            if python:
                results_returned = True
                result_code = result
                result = execute_python(result)
                if (result is not None) and (not result): # Execution complete succesfully, but result was empty
                    results_returned = False
                    result_code = "NOTE: Code returned no results\n\n" + result_code
                    
                    prev_iteration_no_result_notification = f"Note: Task '{task}' completed but returned no results. Please decide if you should retry. If so, please change something so that you will get a result."
                    print(Fore.BLUE + f"Task '{task}' completed but returned no results")

                if "MYGENE" in task:
                    params = get_code_params(
                        result_code,
                        preparam_text="mygene.MyGeneInfo()",
                        postparam_text="gene_results = mg.query(",
                    )
                    if results_returned:
                        cache["MYGENE"].append(f"---\n{params}---\n")
                    else:
                        cache["MYGENE"].append(f"---\nNote: This call returned no results\n{params}---\n")
                    result = process_mygene_result(result)

                if "PUBMED" in task:
                    params = get_code_params(
                        result_code,
                        preparam_text="from Bio import Entrez",
                        postparam_text="search_handle = Entrez.esearch(",
                    )
                    if results_returned:
                        cache["PUBMED"].append(f"---\n{params}---\n")
                    else:
                        cache["PUBMED"].append(f"---\nNote: This call returned no results\n{params}---\n")
                    result = process_pubmed_result(result)

            if type(result) is not list:
                result = [result]

            doc_store_key = str(task_id_counter) + "_" + task
            doc_store["tasks"][doc_store_key] = {}
            doc_store["tasks"][doc_store_key]["results"] = []

            if result_code:
                doc_store["tasks"][doc_store_key]["result_code"] = result_code

            for i, r in enumerate(result):

                # We have result, clear out no result notification
                prev_iteration_no_result_notification = ""
                r = str(r)[
                    :RESULT_CUTOFF
                ]  # Occasionally an enormous result will slow the program to a halt. Not ideal to lose results but putting in place for now.
                vectorized_data = get_ada_embedding(r)
                task_id = f"doc_id_{task_id_counter}_{i}"
                insert_doc_llama_index(index, vectorized_data, task_id, r)

                doc_store["tasks"][doc_store_key]["results"].append(
                    {
                        "task_id_counter": task_id_counter,
                        "vectorized_data": vectorized_data,
                        "output": r,
                    }
                )

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
MAX_ITERATIONS = 3
OBJECTIVE = "Research possible connections to arythmia and deafness."


# If you would like to reload a previous state, comment out run(OBJECTIVE=OBJECTIVE, MAX_ITERATIONS=MAX_ITERATIONS, TOOLS=TOOLS) and uncomment #run(reload_path="out/Cure breast cancer_2023-04-25_16-38-42")
# Then put your path in to your saved state.
# Note that state reloading is not backwards compatible. Saved states before this change cannot be reloaded.

# Fresh Run
run(OBJECTIVE=OBJECTIVE, MAX_ITERATIONS=MAX_ITERATIONS, TOOLS=TOOLS)

# Reload state and resume run
# TOOLS and MAX_ITERATIONS can also be passed in when reloading state.
# run(reload_path="out/Cure breast cancer_2023-04-25_16-38-42")
