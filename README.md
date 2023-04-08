# INSIGHT

## How to run

Sign up for pinecone (free) and OpenAI (also free but will cost a few cents to use)

Expose the following environment vars
    - EMAIL
    - PINECONE_API_KEY
    - PINECONE_ENV
    - OPENAI_API_KEY
    - OPENAI_ORG

run `pip install -r requirements.txt`
run `python main.py`

BE MINDFUL OF EXPENSES!!
    - Currently an execution for a few minutes should cost no more than a few cents. This will go up if you use a more powerful model like GPT-4

NOTE:
    At the bottom of the main.py loop there is a break statement to safeguard against the loop running forever.
        ```
        if task_id_counter > 15:
            break
        ```
    The program will cease execution after 15 iterations and your state will be lost. This is fine while we continue developement, but once the program starts to produce useful outputs this should change.


## TODOs

Explicitly inform boss about worker errors. Ensure boss handles them gracefully.
    - example: code execution issue
Periodically consolidate memories
    - Maybe llama index will let us do something like this out of the box. In which case we wouldn't even need pinecone
Expand APIs
    - more useful examples
    - different APIs
    - NOTE: prompt space is out most scarce commodity. Be as terse as possible and be careful not to add redundancies
Implement token limit constraints
    - currently it's possible to pass an agent a prompt with say 4096 tokens, leaving it space for only one single token for its completion. This should be fixed, but whatever error is passed should be propagated and fixed gracefully
Rework project structure
    - potentially abstractions as well?
Allow way to save and reload program
The pinecone state persists across executions. This should only happen if we explicitly save the state. Llama index and the gpt sessions do not persist.
Implement logging
    - Ideally the output of every run should be saved to its own log file