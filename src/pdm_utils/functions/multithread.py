"""
Functions to multithread process a list of inputs.
"""

from queue import Queue
import threading


class MixedThread(threading.Thread):
    def __init__(self, thread_id, target, work_queue, queue_lock):
        threading.Thread.__init__(self)

        self.name = thread_id

        # Function for thread to call
        self.target = target
        # Queue filled with tuples containing target function args
        self.queue = work_queue
        # Lock controlling access to queue accross multiple threads
        self.lock = queue_lock

    def run(self):
        while not self.stack.empty():
            # Tell other threads "I'm accessing the queue right now"
            self.lock.acquire()

            work_item = self.queue.get()
            # Let other threads have a turn
            self.lock.release()

            # Run target function with args retrieved from queue
            self.target(*work_item)


def create_threads(target, work_stack, stack_lock, num_threads):
    """Function to create a pseudo-MixIn Thread class for multithreading.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :param work_stack: Stack containing tuples of target function args.
    :type work_stack: LifoQueue
    :param stack_lock: Lock for controlling access of stack accross threads.
    :type stack_lock: Lock
    :param num_threads: Number of threads to be created to process work stack.
    :type num_threads: int
    :returns: Returns a list of MixIn Thread objects with target function
    :rtype: list[Thread]
    """
    threads = []
    for x in range(num_threads):
        threads.append(MixedThread(x+1))

    return threads


def multithread(target, work_items, threads=1):
    """Runs list of work items with specified number of threads.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :type work_items: List containing tuples of target function args.
    :type work_stack: list
    :param threads: Number of threads to be created to process work_stack.
    :type threads: int
    """
    # Create lock which prevents queue access deadlocks
    lock = threading.Lock()
    queue = Queue()

    for item in work_items:
        queue.put(item)

    threads = create_threads(target, queue, lock, threads)

    for thread in threads:
        # Calls thread run()
        thread.start()

    for thread in threads:
        # Blocks and waits for threads to finish
        thread.join()
