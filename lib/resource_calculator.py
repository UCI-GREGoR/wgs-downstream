import random


def select_queue(queuename: str, all_queues: dict) -> str:
    """
    use a yaml string representation of a queue group and a configured
    yaml array of queue names to select queue targets for a rule resource
    """
    if queuename in all_queues.keys():
        return random.choice(all_queues[queuename])
    raise ValueError(
        "Configured queue set does not match anything in user resource config: {}".format(
            queuename
        )
    )
