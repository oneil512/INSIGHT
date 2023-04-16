import os

from colorama import Fore
from collections import deque, defaultdict
from Bio import Entrez
import xml.etree.ElementTree as ET

import logging
import llama_index
from llama_index import GPTSimpleVectorIndex


from agents import boss_agent, worker_agent
from utils import execute_python, get_ada_embedding, insert_doc_llama_index, get_code_params, query_knowledge_base

logging.getLogger('llama_index').setLevel(logging.WARNING)

Entrez.email = os.environ['EMAIL']

MAX_TOKENS = 4097
OBJECTIVE = "Cure breast cancer"
tools = ["MYGENE", "PUBMED"]

# Create llama index
llama_index = GPTSimpleVectorIndex([])

task_id_counter = 1
task_list = deque()

tool_description = """
1) Query mygene API. This is useful for finding information on specific genes, or genes associated with the search query. If you wish to make a task to create an API request to mygene then simply say 'MYGENE:' followed by what you would like to search for. Example: 'MYGENE: look up information on genes that are linked to cancer'
2) Query PubMed API. This is useful for searching biomedical literature and studies on any medical subject. If you wish to make a task to create an API request to the PubMed API then simply say 'PUBMED:' followed by what you would like to search for. Example: 'PUBMED: Find recent developments in HIV research'
""".strip()

completed_tasks = []

print(Fore.CYAN, '\n*****OBJECTIVE*****\n')
print(OBJECTIVE)

cache = defaultdict(list)

while True:
    result_code = None
    python = False
    cache_params = None

    if task_id_counter > 1:
        executive_summary = query_knowledge_base(llama_index)
    else:
        executive_summary = "No tasks completed yet."

    task_list, thoughts = boss_agent(
        objective=OBJECTIVE,
        tool_description=tool_description,
        task_list=task_list,
        executive_summary=executive_summary,
        completed_tasks=completed_tasks
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
        if task_id_counter > 1:
            context = query_knowledge_base(llama_index, query=f"Provide as much useful context as possible for this task: {task}")

        if any(tool in task for tool in tools):
            python = True

        if 'MYGENE' in task and 'MYGENE' in cache:
            cache_params = cache['MYGENE']
        
        if 'PUBMED' in task and 'PUBMED' in cache:
            cache_params = cache['PUBMED']

        print(Fore.RED + "\n*****NEXT TASK*****\n")
        print("task id: ", task_id_counter, 'task: ', task)

        result = worker_agent(OBJECTIVE, task, context, cache_params, python)
        completed_tasks.append(task)
        print(Fore.GREEN + '\n*****TASK RESULT*****\n')
        print(Fore.GREEN + result)

        if python:
            result_code = result
            result = execute_python(result)

        if 'PUBMED' in task:
            root = ET.fromstring(result)
            result = []
            for AbstractText in root.iter('AbstractText'):
                result.append(AbstractText.text)

        if type(result) is not list:
            result = [result]

        # Store data
        if python:
            if 'MYGENE' in task:
                params = get_code_params(result_code, preparam_text = "mygene.MyGeneInfo()", postparam_text = "mg.query(")
                cache['MYGENE'].append(params)

            if 'PUBMED' in task:
                params = get_code_params(result_code, preparam_text = "from Bio import Entrez", postparam_text = "search_handle = Entrez.esearch(")
                cache['PUBMED'].append(params)

        for i, r in enumerate(result):
            vectorized_data = get_ada_embedding(str(r))
            task_id = f"doc_id_{task_id_counter}_{i}"
            metadata = {"Task": task, "Result": str(r)}

            insert_doc_llama_index(llama_index, vectorized_data, task_id, str(r))


        task_id_counter += 1

    if task_id_counter > 15:
        break