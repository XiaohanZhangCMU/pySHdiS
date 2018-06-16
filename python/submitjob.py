jobName = 'abq2paraInput'
mdb.Job(name=jobName, model='Beam')
#, description='', type=ANALYSIS, atTime=None, 
#    waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, 
#    getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
#    nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
#    contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
#    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

myJob = mdb.Job(name=jobName, model='Beam')
myJob.submit() 
myJob.waitForCompletion()
