import logging
import math
import os
import time
from collections import defaultdict, deque

from Bio import Entrez
from colorama import Fore, Style
from pyfiglet import Figlet

from agents import boss_agent, worker_agent
from config import EMAIL, OPENAI_API_KEY
from interface import prompt_user
from utils import (
    get_key_results,
    insert_doc_llama_index,
    load,
    query_knowledge_base,
    save,
    handle_python_result,
    handle_results,
    create_index,
    read_file,
    select_task,
)

logging.getLogger("llama_index").setLevel(logging.WARNING)

Entrez.email = EMAIL
MAX_TOKENS = 4097
api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]


def run_(
    api_key,
    doc_store,
    OBJECTIVE,
    RESULT_CUTOFF,
    TOOLS,
    master_index,
    task_id_counter,
    completed_tasks,
    cache,
    tool_description,
    previous_result,
    previous_task,
    task_list,
    summaries,
):
    
    index = create_index(api_key)

    task_list = boss_agent(
        objective=OBJECTIVE,
        tool_description=tool_description,
        task_list=task_list,
        summaries=summaries,
        completed_tasks=completed_tasks,
        previous_task=previous_task,
        previous_result=previous_result
    )

    if not task_list:
        print(Fore.RED + "TASK LIST EMPTY")
        return

    print(Fore.WHITE + "\033[1m\n*****TASK LIST*****\n\033[0m")
    for i, t in enumerate(task_list):
        print(f"{i + 1}) {t}")

    task, task_list = select_task(task_list)

    print(Fore.RED + "\033[1m\n*****NEXT TASK*****\n\033[0m")
    print("task: ", task)

    result, result_is_python = worker_agent(OBJECTIVE, task, master_index, cache, TOOLS)
    completed_tasks.append(task)

    print(Fore.GREEN + "\033[1m\n*****TASK RESULT*****\n\033[0m")
    print(Fore.GREEN + result)

    doc_store_task_key = str(task_id_counter) + "_" + task
    doc_store["tasks"][doc_store_task_key] = {}
    doc_store["tasks"][doc_store_task_key]["results"] = []

    if result_is_python:
        result = handle_python_result(result, cache, task, doc_store, doc_store_task_key)

    if result:
        handle_results(result, index, doc_store, doc_store_task_key, task_id_counter, RESULT_CUTOFF)

        if index.docstore.docs:
            executive_summary, citation_data = query_knowledge_base(index, list_index=False)
            insert_doc_llama_index(index=master_index, doc_id=str(task_id_counter), data=executive_summary, metadata={"citation_data": citation_data})
            doc_store["tasks"][doc_store_task_key]["executive_summary"] = executive_summary
            summaries.append(executive_summary)
            index = create_index(api_key=api_key)

    return result, task, task_list, summaries

    


def run(
    api_key,
    OBJECTIVE="",
    RESULT_CUTOFF=20000,
    MAX_ITERATIONS=1,
    TOOLS=["MYGENE", "PUBMED", "MYVARIANT"],
    master_index=None,
    task_id_counter=1,
    task_list=deque(),
    completed_tasks=[],
    cache=defaultdict(list),
    current_datetime="",
    reload_path="",
    my_data_path="",
):
    start_time = time.time()
    reload_count = 0
    summaries = []
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
                summaries,
            ) = load(reload_path)
            MAX_ITERATIONS += task_id_counter - 1
        except Exception as e:
            print(f"Issue loading state {e}, with path {reload_path}")
            raise Exception("Cannot reload state.")

    else:
        if not master_index:
            master_index = create_index(api_key=api_key)
        current_datetime = current_datetime or str(time.strftime("%Y-%m-%d_%H-%M-%S"))

    if my_data_path:
        my_data = read_file(my_data_path)

        temp_index = create_index(api_key=api_key)
        insert_doc_llama_index(temp_index, data=my_data, doc_id="my_data")
        executive_summary, _ = query_knowledge_base(temp_index, list_index=False)

        insert_doc_llama_index(index=master_index, doc_id="my_data", data=executive_summary)
        summaries.append(executive_summary)

    doc_store = {"tasks": {}}

    tool_description_mapping = {
        "MYGENE": """1) Query mygene API. This is useful for finding information on specific genes, or genes associated with the search query. If you wish to make a task to create an API request to mygene then simply say 'MYGENE:' followed by what you would like to search for. Example: 'MYGENE: look up information on genes that are linked to cancer'""",
        "PUBMED": """2) Query PubMed API. This is useful for searching biomedical literature and studies on any medical subject. If you wish to make a task to create an API request to the PubMed API then simply say 'PUBMED:' followed by what you would like to search for. Example: 'PUBMED: Find recent developments in HIV research'""",
        "MYVARIANT": """3) Query myvariant API. This is useful for finding information on specific genetic variants. If you wish to make a task to create an API request to myvariant then simply say 'MYVARIANT:' followed by the specific genetic variant you are interested in. You can specify by rsID, ClinVar, or in a standardized format for describing a genetic variant like 'chr1:g.35367G>A'. Example: 'MYVARIANT: look up information on the variant chr1:g.35367G>A'""",
    }

    tool_description = ""

    for tool in TOOLS:
        tool_description += f"{tool_description_mapping[tool]}\n"

    print(Fore.CYAN, "\033[1m\n*****OBJECTIVE*****\n\033[0m")
    print(OBJECTIVE)


    result, completed_tasks = [], []
    task = ""
    task_list = deque()
    cache=defaultdict(list)

    for _ in range(MAX_ITERATIONS):
        result, task, task_list, summaries = run_(
            api_key=api_key,
            doc_store=doc_store,
            OBJECTIVE=OBJECTIVE,
            RESULT_CUTOFF=RESULT_CUTOFF,
            TOOLS=TOOLS,
            master_index=master_index,
            task_id_counter=task_id_counter,
            task_list=task_list,
            completed_tasks=completed_tasks,
            cache=cache,
            previous_task=task,
            tool_description=tool_description,
            previous_result=result,
            summaries=summaries,
        )
        
        task_id_counter += 1
        

    doc_store["key_results"] = get_key_results(master_index, OBJECTIVE, top_k=20)

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
        summaries
    )

    end_time = time.time()
    total_time = end_time - start_time
    print(Fore.RED + f"Total run time: {total_time:.2f} seconds")


if __name__ == "__main__":
    fig = Figlet(font='slant')
    title_art = fig.renderText("INSIGHT")
    print(Fore.LIGHTMAGENTA_EX + f"\033[1m{title_art}\033[0m")

    objective, tool_flags, iterations, reload_path, my_data_path = prompt_user()
    tools = [key for key, value in tool_flags.items() if value]

    run(api_key=api_key, OBJECTIVE=objective, MAX_ITERATIONS=iterations, TOOLS=tools, my_data_path=my_data_path, reload_path=reload_path)