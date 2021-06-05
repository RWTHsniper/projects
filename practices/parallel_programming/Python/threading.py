# -*- coding: utf-8 -*-
"""
Created on Thu May  6 20:54:57 2021

@author: golde
"""

import threading
import time


class Worker(threading.Thread):
    def __init__(self, name):
        super().__init__()
        self.name = name            # thread 이름 지정

    def run(self):
        print("sub thread start ", threading.currentThread().getName())
        time.sleep(3)
        print("sub thread end ", threading.currentThread().getName())


print("main thread start")
for i in range(5):
    name = "thread {}".format(i)
    t = Worker(name)                # sub thread 생성
    t.start()                       # sub thread의 run 메서드를 호출

print("main thread end")



# 반복문을 통해 여러 서브 스레드를 생성해야하는 경우에는 생성된 스레드 객체를 파이썬 리스트에 저장한 후 반복문을 이용해서 각 객체에서 join( ) 메서드를 호출할 수 있습니다.


class Worker(threading.Thread):
    def __init__(self, name):
        super().__init__()
        self.name = name            # thread 이름 지정

    def run(self):
        print("sub thread start ", threading.currentThread().getName())
        time.sleep(5)
        print("sub thread end ", threading.currentThread().getName())


print("main thread start")

threads = []
for i in range(3):
    thread = Worker(i)
    thread.start()              # sub thread의 run 메서드를 호출
    threads.append(thread)


for thread in threads:
    thread.join()

print("main thread post job")
print("main thread end")