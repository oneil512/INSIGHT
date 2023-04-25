import logging
import os
from functools import partial

import backoff
import llama_index
import json
import markdown
import openai
import tiktoken
from colorama import Fore
from llama_index import Document
from llama_index.indices.composability import ComposableGraph
import xml.etree.ElementTree as ET

from api.mygene_api import mygene_api
from api.pubmed_api import pubmed_api
from config import OPENAI_API_KEY

logging.getLogger("llama_index").setLevel(logging.WARNING)


MAX_TOKENS = 4097

api_info_mapping = {"mygene": mygene_api, "PubMed": pubmed_api}

openai.api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]


def num_tokens_from_string(string: str, encoding_name: str = "gpt2") -> int:
    """Returns the number of tokens in a text string."""

    encoding = tiktoken.get_encoding(encoding_name)
    num_tokens = len(encoding.encode(string))
    return num_tokens


def get_key_results(index, objective, top_k=50):
    """Run final queries over retrieved documents and store in doc_store."""

    print(Fore.CYAN + "\n*****COMPILING KEY RESULTS*****\n")

    key_results = []

    queries = [
        "Give a brief high level summary of all the data. Cite your sources with the citation information.",
        "Briefly list all the main points that the data covers. Cite your sources with the citation information.",
        "Give all of the key insights about the data. Cite your sources with the citation information.",
        "Generate several creative hypotheses given the data.",
        "What are some high level research directions to explore further given the data?",
        "Describe the key findings in great detail. Do not include filler words. Cite your sources with the citation information.",
        f"Do your best to answer the objective: {objective} given the information. Cite your sources with the citation information."
    ]

    for query in queries:
        print(Fore.CYAN + f"\nCOMPILING RESULT {query}\n")
        res = None
        try:
            res = query_knowledge_base(index=index, query=query, top_k=top_k)
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

    for res in result:

        json_data = res

        _id = json_data.get('_id')
        _version = json_data.get('_version')
        name = json_data.get('name')
        refseq = json_data.get('refseq', {}).get('genomic', [])
        symbol = json_data.get('symbol')
        taxid = json_data.get('taxid')
        pathway = json_data.get('pathway')
        type_of_gene = json_data.get('type_of_gene')
        summary = json_data.get('summary')
        kegg = json_data.get('kegg', [])
        pid = json_data.get('pid', {})
        reactome = json_data.get('reactome', [])
        wikipathways = json_data.get('wikipathways', {})

        output = f"ID: {_id}\n"
        output += f"\nVersion: {_version}\n"
        if name:
            output += f"Name: {name}\n"
        if refseq:
            output += f"RefSeq: {', '.join(refseq)}\n"
        if symbol:
            output += f"Symbol: {symbol}\n"
        if taxid:
            output += f"Tax ID: {taxid}\n"

        output += f"PATHWAYS\n\n"
        output += "Kegg:\n"
        for item in kegg:
            output += f" ID: {item.get('id', '')}"
            output += f" Name: {item.get('name', '')}"

        output += "\nPid:\n"
        output += f" ID: {pid.get('id', '')}"
        output += f" Name: {pid.get('name', '')}"

        output += "\nReactome:\n"
        for item in reactome:
            output += f" ID: {item.get('id', '')}"
            output += f" Name: {item.get('name', '')}"

        output += "\nWikipathways:\n"
        output += f"  ID: {wikipathways.get('id', '')}"
        output += f"  Name: {wikipathways.get('name', '')}"
        if type_of_gene:
            output += f"Type of gene: {type_of_gene}\n"
        if summary:
            output += f"Summary of {name}: {summary}\n"
        output += '\n'

        processed_result.append(output)

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

        processed_result.append(res_)

    return processed_result

def prune_gene_results(res):
    pass


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

The example doesn't have to be followed exactly. You should change it to fit your specific task.

        """.strip()

    return prompt


def get_ada_embedding(text):
    ada_embedding_max_size = 8191
    text = text.replace("\n", " ")

    if len(text) > ada_embedding_max_size:
        # There must be a better way to do this.
        text = text[:ada_embedding_max_size]
    return openai.Embedding.create(input=[text], model="text-embedding-ada-002")[
        "data"
    ][0]["embedding"]


def insert_doc_llama_index(index, embedding, doc_id, metadata):
    doc = Document(text=metadata, embedding=embedding, doc_id=doc_id)
    index.insert(doc)


def query_knowledge_base(
    index,
    query="Give a detailed but terse overview of all the information. Start with a high level summary and then go into details. Do not include any further instruction. Do not include filler words. Do not include citation information.",
    response_mode="tree_summarize",
    top_k=50,
):
    # From llama index docs: Empirically, setting response_mode="tree_summarize" also leads to better summarization results.
    query_response = index.query(
        query, similarity_top_k=top_k, response_mode=response_mode
    )
    return query_response.response


def parser(instruction, content):
    return (
        openai.Completion.create(
            engine="text-davinci-003",
            prompt=instruction + "\nHere is the content to parse:\n" + content,
            temperature=0.0,
        )
        .choices[0]
        .text.strip()
    )


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


def save(index, doc_store, OBJECTIVE, current_datetime):
    path = os.path.join("./out", OBJECTIVE + "_" + current_datetime)
    make_dir(path)

    index.save_to_disk(os.path.join(path, 'index.json'))

    if "key_results" in doc_store:
        header = f"# {OBJECTIVE}\nDate: {current_datetime}\n\n"
        write_file(os.path.join(path, "key_findings.md"), header, mode="a+")
        for res in doc_store["key_results"]:
            content = f"{res[0]}{res[1]}"
            write_file(os.path.join(path, "key_findings.md"), content, mode="a+")

    for task, doc in doc_store["tasks"].items():
        doc_path = os.path.join(path, task)
        make_dir(doc_path)
        result_path = os.path.join(doc_path, "results")
        make_dir(result_path)

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

def load():
    pass
