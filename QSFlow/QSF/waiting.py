import time


class waiting:
    def __init__(self, name,dir):
        self.dir = dir
        self.name = name
        made = False
        while made == False:
            with open(
                    f'{self.dir}/{self.name}script/{self.name}-final.log',
                    'r') as outputF:
                a = outputF.readlines()

                if len(a) == 0:
                    i =0
                    while len(a) ==0:
                        a= outputF.readlines()
                        i = i+1
                        print("log file is courpted, retrying")
                        if i ==10:
                            print("please retry, terminating program")
                            break
                line = a[-1].strip()
                if line[0:33] == "Normal termination of Gaussian 16":
                    made = True
                    print("Now trasfering charges")
                else:
                    print("DFT not done, wating for 10 sec befroe rechecking")
                    time.sleep(10)
