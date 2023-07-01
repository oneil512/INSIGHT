import os
from config import EMAIL, OPENAI_API_KEY
from utils import query_knowledge_base
import time
from langchain import OpenAI
from llama_index import GPTVectorStoreIndex, LLMPredictor, ServiceContext


def query_data(index):
    while True:
        try:
            query = input("Please enter your question (or 'exit' to quit): ").strip()
            
            # If the user wants to exit, break the loop
            if query.lower() == "exit":
                break

            reponse, citation_data = query_knowledge_base(
                index,
                query=query,
                response_mode="tree_summarize",
                top_k=50,
                list_index=False
            )

            print(reponse, citation_data)

        except Exception as e:
            # If there is any error in the user's input, print an error message
            # but keep the loop going
            print(f"An error occurred: {str(e)}. Please try again.")

        # To avoid potential endless loops, a small delay is introduced
        time.sleep(0.1)

if __name__ == "__main__":
    api_key = OPENAI_API_KEY or os.environ["OPENAI_API_KEY"]

    PATH_TO_INDEX = 'out/Cure breast cancer_2023-04-29_14-58-41/index.json'

    llm_predictor = LLMPredictor(
        llm=OpenAI(
            temperature=0,
            openai_api_key=api_key,
            model_name="text-davinci-003",
            max_tokens=2000,
        )
    )
    service_context = ServiceContext.from_defaults(llm_predictor=llm_predictor)
    index = GPTVectorStoreIndex.load_from_disk(
        PATH_TO_INDEX, service_context=service_context
    )

    query_data(index)