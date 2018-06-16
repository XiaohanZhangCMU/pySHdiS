import time
import os.path
while os.path.exists('instructions.py'):
    execfile('instructions.py')
    time.sleep(1)
