# INSIGHT

## How to run

Sign up for pinecone (free) and OpenAI (also free but will cost a few cents to use)

Expose the following environment vars<br>
    - EMAIL<br>
    - PINECONE_API_KEY<br>
    - PINECONE_ENV<br>
    - OPENAI_API_KEY<br>
    - OPENAI_ORG<br>

run `pip install -r requirements.txt`<br>
run `python main.py`

BE MINDFUL OF EXPENSES!!<br>
    - Currently an execution for a few minutes should cost no more than a few cents. This will go up if you use a more powerful model like GPT-4<br>

NOTE:<br>
    At the bottom of the main.py loop there is a break statement to safeguard against the loop running forever.<br>
        ```py
        if task_id_counter > 15:
            break
        ```
    The program will cease execution after 15 iterations and your state will be lost. This is fine while we continue developement, but once the program starts to produce useful outputs this should change.<br>


## TODOs

Explicitly inform boss about worker errors. Ensure boss handles them gracefully.<br>
    - example: code execution issue<br><br>
Periodically consolidate memories<br>
    - Maybe llama index will let us do something like this out of the box. In which case we wouldn't even need pinecone<br><br>
Expand APIs<br>
    - more useful examples<br>
    - different APIs<br>
    - NOTE: prompt space is out most scarce commodity. Be as terse as possible and be careful not to add redundancies<br><br>
Implement token limit constraints<br>
    - currently it's possible to pass an agent a prompt with say 4096 tokens, leaving it space for only one single token for its completion. This should be fixed, but whatever error is passed should be propagated and fixed gracefully<br><br>
Rework project structure<br>
    - potentially abstractions as well?<br>
Allow way to save and reload program<br>
The pinecone state persists across executions. This should only happen if we explicitly save the state. Llama index and the gpt sessions do not persist.<br><br>
Implement logging<br>
    - Ideally the output of every run should be saved to its own log file<br>