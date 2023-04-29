import json
import logging
import os
import time
import xml.etree.ElementTree as ET
from collections import defaultdict, deque
from functools import partial

import backoff
import llama_index
import markdown
import openai
import tiktoken
from colorama import Fore
from langchain import OpenAI
from llama_index import Document, GPTSimpleVectorIndex, LLMPredictor, ServiceContext, GPTListIndex
from llama_index.indices.composability import ComposableGraph

from api.mygene_api import mygene_api
from api.pubmed_api import pubmed_api
from config import OPENAI_API_KEY

logging.getLogger("llama_index").setLevel(logging.WARNING)

#file_handler = logging.FileHandler('utils.log')
# Configure the logging settings
#logging.basicConfig(level=logging.INFO, handlers=[file_handler])


MAX_TOKENS = 4097
api_info_mapping = {"mygene": mygene_api, "PubMed": pubmed_api}

api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]
openai.api_key = api_key


def num_tokens_from_string(string: str, encoding_name: str = "gpt2") -> int:
    """Returns the number of tokens in a text string."""

    encoding = tiktoken.get_encoding(encoding_name)
    num_tokens = len(encoding.encode(string))
    return num_tokens


def get_key_results(index, objective, top_k=20):
    """Run final queries over retrieved documents and store in doc_store."""

    if not index.docstore.docs:
        print(Fore.RED + "\n! NO TASKS RETURNED RESULTS. PLEASE TWEAK YOUR OBJECTIVE AND CHECK SPELLING. !\n")
        return []


    print(Fore.CYAN + "\n*****COMPILING KEY RESULTS*****\n")

    key_results = []

    queries = [
        "Give a brief high level summary of all the data. Cite your sources with the citation information. Do not make up citation information.",
        "Briefly list all the main points that the data covers. Cite your sources with the citation information. Do not make up citation information.",
        "Give all of the key insights about the data. Cite your sources with the citation information. Do not make up citation information.",
        "Generate several creative hypotheses given the data.",
        "What are some high level research directions to explore further given the data?",
        "Describe the key findings in great detail. Do not include filler words. Cite your sources with the citation information. Do not make up citation information.",
        f"Do your best to answer the objective: {objective} given the information. Cite your sources with the citation information. Do not make up citation information.",
    ]

    for query in queries:
        print(Fore.CYAN + f"\nCOMPILING RESULT {query}\n")
        res = None
        try:
            res = query_knowledge_base(index=index, query=query, list_index=False, top_k=top_k)
        except Exception as e:
            print(f"Exception getting key result {query}, error {e}")

        if res:
            query = f"## {query}\n\n"
            html = markdown.markdown(res)
            key_results.append((query, f"{html}\n\n\n\n"))

    print(Fore.CYAN + f"\nRESULTS COMPILED. SAVED TO DIRECTORY `out`\n")

    return key_results


def get_max_completion_len(prompt):
    tokens = num_tokens_from_string(prompt)
    return MAX_TOKENS - tokens


def execute_python(code: str):
    # ret is defined in the code string
    loc = {}
    try:
        exec(code, globals(), loc)

    except Exception as e:
        print(f"Exception executing code {code}, {e}")
        return

    return loc["ret"]


def process_mygene_result(result):
    processed_result = []

    # Each result will be split into 2 documents: summary and pathway
    for json_data in result:

        name = json_data.get("name")
        refseq_genomic = json_data.get("refseq", {}).get("genomic", [])
        refseq_rna = json_data.get("refseq", {}).get("rna", [])
        symbol = json_data.get("symbol")
        taxid = json_data.get("taxid")
        type_of_gene = json_data.get("type_of_gene")
        pos = json_data.get("genomic_pos_hg19")
        summary = json_data.get("summary")

        output_summary = ""
        # Summary
        if name:
            output_summary += f"Gene Name: {name}\n"
        if refseq_genomic:
            output_summary += f"RefSeq genomic: {', '.join(refseq_genomic)}\n"
        if refseq_rna:
            output_summary += f"RefSeq rna: {', '.join(refseq_rna)}\n"
        if symbol:
            output_summary += f"Symbol: {symbol}\n"
        if taxid:
            output_summary += f"Tax ID: {taxid}\n"
        if type_of_gene and type_of_gene != 'unknown':
            output_summary += f"Type of gene: {type_of_gene}\n"
        if pos:
            output_summary += f"Position: {pos}\n"
        if summary:
            output_summary += f"Summary of {name}: {summary}\n"

        output_summary = output_summary.strip()

        #logging.info(f"Mygene Summary result {name}, length is {str(len(output_summary))}")
        if output_summary:
            processed_result.append(output_summary)
        

        # Pathway
        pathway = json_data.get("pathway")
        if pathway:
            kegg = pathway.get("kegg", [])
            pid = pathway.get("pid", [])
            reactome = pathway.get("reactome", [])
            wikipathways = pathway.get("wikipathways", [])
            netpath = pathway.get("netpath", [])
            biocarta = pathway.get("biocarta", [])

            pathway_elements = {"kegg": kegg, "pid": pid, "reactome": reactome, "wikipathways": wikipathways, "netpath": netpath, "biocarta": biocarta}

            # mygene returns dicts instead of lists if singleton
            # Wrap with list if not list
            for k,v in pathway_elements.items():
                if type(v) is not list:
                    pathway_elements[k] = [v]


            output_pathway = ""
            if name:
                output_pathway += f"Gene Name: {name}\n"
            if symbol:
                output_pathway += f"Symbol: {symbol}\n"
            if taxid:
                output_pathway += f"Tax ID: {taxid}\n"
            if type_of_gene and type_of_gene != 'unknown':
                output_pathway += f"Type of gene: {type_of_gene}\n"
            if refseq_genomic:
                output_pathway += f"RefSeq genomic: {', '.join(refseq_genomic)}\n"
            if refseq_rna:
                output_pathway += f"RefSeq rna: {', '.join(refseq_rna)}\n"
            if pos:
                output_pathway += f"Position: {pos}\n"

            output_pathway += f"PATHWAYS\n\n"

            for k,v in pathway_elements.items():
                output_pathway += f"\n{k}:\n"
                for item in v:
                    output_pathway += f" ID: {item.get('id', '')}"
                    output_pathway += f" Name: {item.get('name', '')}"

            #logging.info(f"Mygene Pathway result {name}, length is {len(output_pathway)}")

            output_pathway = output_pathway.strip()
            if output_pathway:
                processed_result.append(output_pathway)

    return processed_result


def process_pubmed_result(result):
    root = ET.fromstring(result)
    processed_result = []

    for article in root:
        res_ = ""
        for title in article.iter("Title"):
            res_ += f"{title.text}\n"
        for abstract in article.iter("AbstractText"):
            res_ += f"{abstract.text}\n"
        for author in article.iter("Author"):
            try:
                res_ += f"{author.find('LastName').text}"
                res_ += f", {author.find('ForeName').text}\n"
            except:
                pass
        for journal in article.iter("Journal"):
            res_ += f"{journal.find('Title').text}\n"
        for volume in article.iter("Volume"):
            res_ += f"{volume.text}\n"
        for issue in article.iter("Issue"):
            res_ += f"{issue.text}\n"
        for pubdate in article.iter("PubDate"):
            try:
                year = pubdate.find("Year").text
                res_ += f"{year}"
                month = pubdate.find("Month").text
                res_ += f"-{month}"
                day = pubdate.find("Day").text
                res_ += f"-{day}\n"
            except:
                pass
        for doi in article.iter("ELocationID"):
            if doi.get("EIdType") == "doi":
                res_ += f"{doi.text}\n"

        if res_:
            processed_result.append(res_)

    return processed_result


def get_code_params(code: str, preparam_text: str, postparam_text: str):
    l = len(preparam_text)

    preparam_index = code.find(preparam_text)
    postparam_index = code.find(postparam_text)

    if preparam_index == -1 or postparam_index == -1:
        return

    params = code[preparam_index + l : postparam_index].strip()

    if params == "":
        return

    return params


def validate_llm_response(goal, response):
    validation_prompt = f"I gave an LLM this goal: '{goal}' and it gave this response: '{response}'. Is this reasonable, or did something go wrong? [yes|no]"
    validation_response = (
        openai.Completion.create(
            engine="text-davinci-003", prompt=validation_prompt, temperature=0.0
        )
        .choices[0]
        .text.strip()
    )

    if validation_response.lower() == "yes":
        return True
    else:
        return False


def generate_tool_prompt(task):
    if "PubChem" in task:
        api_name = "PubChem"
    elif "MYGENE" in task:
        api_name = "mygene"
    elif "PUBMED" in task:
        api_name = "PubMed"
    else:
        print(f"Error. Tool not found in task: {task}")
        return None

    api_info = api_info_mapping[api_name]

    prompt = f"""You have access to query the {api_name} API. If a task starts with '{api_name.upper()}:' then you should create the code to query the {api_name} API based off the documentation and return the code to complete your task. If you use the {api_name} API, do not answer with words, simply answer with the code to query the API and then cease output. Be sure that it is a valid API call that will execute in a python interpreter.
---
Here is the {api_name} documentation
{api_info}
---

You should change the parameters to fit your specific task.

        """.strip()

    return prompt


def get_ada_embedding(text):
    ada_embedding_max_size = 8191
    text = text.replace("\n", " ")

    if num_tokens_from_string(text) > ada_embedding_max_size:
        # There must be a better way to do this.
        text = text[:ada_embedding_max_size]
    return openai.Embedding.create(input=[text], model="text-embedding-ada-002")[
        "data"
    ][0]["embedding"]


def insert_doc_llama_index(index, doc_id, metadata, embedding=None):
    if not embedding:
        embedding = get_ada_embedding(metadata)
    doc = Document(text=metadata, embedding=embedding, doc_id=doc_id)
    index.insert(doc)


def handle_python_result(result, cache, task, doc_store, doc_store_task_key):

    results_returned = True
    doc_store["tasks"][doc_store_task_key]["result_code"] = result

    executed_result = execute_python(result)

    if type(executed_result) is list:
        executed_result = list(filter(lambda x : x, executed_result))

    if executed_result is None:
        result = "NOTE: Code did not run succesfully\n\n" + result
        print(Fore.BLUE + f"Task '{task}' failed. Code {result} did not run succesfully.")
        if "MYGENE" in task:
            cache["MYGENE"].append(f"---\nNote: This call did not run succesfully\n{params}---\n")
        if "PUBMED" in task:
            cache["PUBMED"].append(f"---\nNote: This call did not run succesfully\n{params}---\n")

        return


    if (executed_result is not None) and (not executed_result): # Execution complete succesfully, but executed result was empty list
        results_returned = False
        result = "NOTE: Code returned no results\n\n" + result
        
        print(Fore.BLUE + f"\nTask '{task}' completed but returned no results")


    if "MYGENE" in task:
        params = get_code_params(
            result,
            preparam_text="mygene.MyGeneInfo()",
            postparam_text="gene_results = mg.query(",
        )
        if results_returned:
            cache["MYGENE"].append(f"---\n{params}---\n")
        else:
            cache["MYGENE"].append(f"---\nNote: This call returned no results\n{params}---\n")
        executed_result = process_mygene_result(executed_result)

    if "PUBMED" in task:
        params = get_code_params(
            result,
            preparam_text="from Bio import Entrez",
            postparam_text="search_handle = Entrez.esearch(",
        )
        if results_returned:
            cache["PUBMED"].append(f"---\n{params}---\n")
        else:
            cache["PUBMED"].append(f"---\nNote: This call returned no results\n{params}---\n")
        executed_result = process_pubmed_result(executed_result)


    return executed_result


def handle_results(result, index, doc_store, doc_store_key, task_id_counter, RESULT_CUTOFF):
    if type(result) is not list:
        result = [result]

    for i, r in enumerate(result):
        r = str(r)[
            :RESULT_CUTOFF
        ]  # Occasionally an enormous result will slow the program to a halt. Not ideal to lose results but putting in place for now.
        vectorized_data = get_ada_embedding(r)
        task_id = f"doc_id_{task_id_counter}_{i}"
        insert_doc_llama_index(index=index, doc_id=task_id, metadata=r, embedding=vectorized_data)

        doc_store["tasks"][doc_store_key]["results"].append(
            {
                "task_id_counter": task_id_counter,
                "vectorized_data": vectorized_data,
                "output": r,
            }
        )


def query_knowledge_base(
    index,
    query="Give a detailed but terse overview of all the information. Start with a high level summary and then go into details. Do not include any further instruction. Do not include filler words. Include citation information if it is present. Do not make up citation information.",
    response_mode="tree_summarize",
    top_k=50,
    list_index=False
):
    if not index.docstore.docs:
        print(Fore.RED + "NO INFORMATION IN LLAMA INDEX")
        return
    
    if list_index:
        query_response = index.query(
            query, response_mode="default"
        )
    else:
        # From llama index docs: Empirically, setting response_mode="tree_summarize" also leads to better summarization results.
        query_response = index.query(
            query, similarity_top_k=top_k, response_mode=response_mode
        )
    return query_response.response


def create_index(api_key,summaries=[], temperature=0.0, model_name="text-davinci-003", max_tokens=2000):
    llm_predictor = LLMPredictor(
        llm=OpenAI(
            temperature=temperature,
            openai_api_key=api_key,
            model_name=model_name,
            max_tokens=max_tokens,
            )
    )
    documents = []
    for i, summary in enumerate(summaries):
        documents.append(Document(text=summary, doc_id=str(i)))

    service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)
    return GPTSimpleVectorIndex(documents, service_context=service_context)


def create_graph_index(api_key, indicies=[], summaries=[], temperature=0.0, model_name="text-davinci-003", max_tokens=2000):
    llm_predictor = LLMPredictor(
        llm=OpenAI(
            temperature=temperature,
            openai_api_key=api_key,
            model_name=model_name,
            max_tokens=max_tokens,
            )
    )
    service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)


    graph = ComposableGraph.from_indices(
        GPTListIndex,
        indicies,
        index_summaries=summaries,
        service_context=service_context
    )

    return graph


def create_list_index(api_key, summaries=[], temperature=0.0, model_name="text-davinci-003", max_tokens=2000):
    llm_predictor = LLMPredictor(
        llm=OpenAI(
            temperature=temperature,
            openai_api_key=api_key,
            model_name=model_name,
            max_tokens=max_tokens,
            )
    )
    service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)
    documents = []
    for i, summary in enumerate(summaries):
        documents.append(Document(text=summary, doc_id=str(i)))
        
    index = GPTListIndex.from_documents(documents, service_context=service_context)
    return index


@backoff.on_exception(
    partial(backoff.expo, max_value=50),
    (openai.error.RateLimitError, openai.error.APIError),
)
def get_gpt_completion(
    prompt,
    temp=0.0,
    engine="text-davinci-003",
    top_p=1,
    frequency_penalty=0,
    presence_penalty=0,
):
    response = openai.Completion.create(
        engine=engine,
        prompt=prompt,
        temperature=temp,
        max_tokens=get_max_completion_len(prompt),
        top_p=top_p,
        frequency_penalty=frequency_penalty,
        presence_penalty=presence_penalty,
    )
    return response.choices[0].text.strip()


@backoff.on_exception(
    partial(backoff.expo, max_value=50),
    (openai.error.RateLimitError, openai.error.APIError),
)
def get_gpt_chat_completion(
    system_prompt, user_prompt, model="gpt-3.5-turbo", temp=0.0
):
    response = openai.ChatCompletion.create(
        model=model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=temp,
    )
    return response.choices[0]["message"]["content"].strip()


### FILE UTILS ###


def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def write_file(path, contents, mode="w"):
    with open(path, mode) as f:
        f.write(contents)


def read_file(path, mode="r"):
    with open(path, mode) as f:
        contents = f.read()

    if not contents:
        print(f"WARNING: file {path} empty")

    return contents


def save(
    index,
    doc_store,
    OBJECTIVE,
    current_datetime,
    task_id_counter,
    task_list,
    completed_tasks,
    cache,
    reload_count,
    summaries,
):
    # Make basepath.
    path = os.path.join("./out", OBJECTIVE + "_" + current_datetime)
    make_dir(path)

    # Save llama index.
    index.save_to_disk(os.path.join(path, "index.json"))

    # Save program state.
    state = {
        "summaries": summaries,
        "reload_count": reload_count,
        "task_id_counter": task_id_counter,
        "task_list": list(task_list),
        "completed_tasks": completed_tasks,
        "cache": dict(cache),
        "current_datetime": current_datetime,
        "objective": OBJECTIVE,
    }
    with open(os.path.join(path, "state.json"), "w") as outfile:
        json.dump(state, outfile)

    # Save results.
    if "key_results" in doc_store:
        if reload_count:
            new_time = str(time.strftime("%Y-%m-%d_%H-%M-%S"))
            header = f"# {OBJECTIVE}\nDate: {new_time}\n\n"
        else:
            header = f"# {OBJECTIVE}\nDate: {current_datetime}\n\n"
        key_findings_path = os.path.join(path, f"key_findings_{reload_count}.md")
        write_file(key_findings_path, header, mode="a+")
        for res in doc_store["key_results"]:
            content = f"{res[0]}{res[1]}"
            write_file(key_findings_path, content, mode="a+")

    for task, doc in doc_store["tasks"].items():
            

        doc_path = os.path.join(path, task)
        make_dir(doc_path)
        result_path = os.path.join(doc_path, "results")
        make_dir(result_path)

        if "executive_summary" in doc:
            write_file(os.path.join(result_path, "executive_summary.txt"), doc["executive_summary"])
        if "result_code" in doc:
            write_file(os.path.join(result_path, "api_call.txt"), doc["result_code"])
        

        for i, result in enumerate(doc["results"]):

            result_path_i = os.path.join(result_path, str(i))
            make_dir(result_path_i)
            write_file(os.path.join(result_path_i, "output.txt"), result["output"])
            write_file(
                os.path.join(result_path_i, "vector.txt"),
                str(result["vectorized_data"]),
            )


def load(path):
    llm_predictor = LLMPredictor(
        llm=OpenAI(
            temperature=0,
            openai_api_key=api_key,
            model_name="text-davinci-003",
            max_tokens=2000,
        )
    )
    service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)
    index = GPTSimpleVectorIndex.load_from_disk(
        os.path.join(path, "index.json"), service_context=service_context
    )
    state_path = os.path.join(path, "state.json")
    if os.path.exists(state_path):
        with open(state_path, "r") as f:
            json_data = json.load(f)

            try:
                reload_count = json_data["reload_count"] + 1
                task_id_counter = json_data["task_id_counter"]
                task_list = json_data["task_list"]
                completed_tasks = json_data["completed_tasks"]
                cache = defaultdict(list, json_data["cache"])
                current_datetime = json_data["current_datetime"]
                objective = json_data["objective"]
                summaries = json_data["summaries"]
            except KeyError as e:
                raise Exception(
                    f"Missing key '{e.args[0]}' in JSON file at path '{state_path}'"
                )

    return (
        index,
        task_id_counter,
        deque(task_list),
        completed_tasks,
        cache,
        current_datetime,
        objective,
        reload_count,
        summaries,
    )
