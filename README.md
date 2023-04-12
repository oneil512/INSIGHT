# INSIGHT

Insight is an autonomous AI that can do medical research. It has a boss agent that takes an objective and some context (to come later) and creates a task list. A worker agent picks up a task from the list and completes it, saving the results to llama index. The boss gets informed of the results and changes/reprioritizes the task list. The workers can call into the pubmed and mygene APIs (more to come). The workers also get context from llama index to help complete their tasks.

Very much a work in progress, but it is showing some early results!

Please reach out to me or contribute if this interests you :)

## How to run

Sign up for OpenAI

Expose the following environment vars<br>
    - EMAIL<br>
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

Right now it seems to get overly focused on one area. Better executive summaries might fix this.
    - Look into composable indicies with llama index


Explicitly inform boss about worker errors. Ensure boss handles them gracefully.<br>
    - example: code execution issue


Expand APIs<br>
    - more useful examples<br>
    - different APIs<br>
    - NOTE: prompt space is out most scarce commodity. Be as terse as possible and be careful not to add redundancies


Implement token limit constraints<br>
    - currently it's possible to pass an agent a prompt with say 4096 tokens, leaving it space for only one single token for its completion. This should be fixed, but whatever error is passed should be propagated and fixed gracefully


Rework project structure<br>
    - potentially abstractions as well?


Allow way to save and reload program


Implement logging<br>
    - Ideally the output of every run should be saved to its own log file<br>