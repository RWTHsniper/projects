# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 23:16:37 2020

@author: golde
"""



from ibapi.client import EClient
from ibapi.wrapper import EWrapper
import threading
import time
import pandas as pd


class TradingApp(EWrapper, EClient):
    def __init__(self):
        EClient.__init__(self,self)
        self.pos_df = pd.DataFrame(columns=['Account', 'Symbol', 'SecType',
                                            'Currency', 'Position', 'Avg cost'])
        self.acc_summary = pd.DataFrame(columns=['ReqId', 'Account', 'Tag', 'Value', 'Currency'])
        self.pnl_summary = pd.DataFrame(columns=['ReqId', 'DailyPnL', 'UnrealizedPnL', 'RealizedPnL'])
        
    def position(self, account, contract, position, avgCost):
        super().position(account, contract, position, avgCost)
        dictionary = {"Account":account, "Symbol": contract.symbol, "SecType": contract.secType,
                      "Currency": contract.currency, "Position": position, "Avg cost": avgCost}
        self.pos_df = self.pos_df.append(dictionary, ignore_index=True)

    def accountSummary(self, reqId, account, tag, value, currency):
        super().accountSummary(reqId, account, tag, value, currency)
        dictionary = {"ReqId":reqId, "Account": account, "Tag": tag, "Value": value, "Currency": currency}
        self.acc_summary = self.acc_summary.append(dictionary, ignore_index=True)
        
    def pnl(self, reqId, dailyPnL, unrealizedPnL, realizedPnL):
        super().pnl(reqId, dailyPnL, unrealizedPnL, realizedPnL)
        dictionary = {"ReqId":reqId, "DailyPnL": dailyPnL, "UnrealizedPnL": unrealizedPnL, "RealizedPnL": realizedPnL}
        self.pnl_summary = self.pnl_summary.append(dictionary, ignore_index=True)
        

def websocket_con():
    app.run()
    
app = TradingApp()      
app.connect("127.0.0.1", 7497, clientId=1)

# starting a separate daemon thread to execute the websocket connection
con_thread = threading.Thread(target=websocket_con, daemon=True)
con_thread.start()
time.sleep(1) # some latency added to ensure that the connection is established


app.reqPositions()
time.sleep(1)
pos_df = app.pos_df
time.sleep(5)

app.reqAccountSummary(1, "All", "$LEDGER:USD") # I want the numbers in USD
time.sleep(1)
acc_summ_df = app.acc_summary

#app.reqPnL(2, "DU111519", "")
app.reqPnL(2, "DU1894924", "")
time.sleep(1)
pnl_summ_df = app.pnl_summary

time.sleep(5)


UnrealizedPnL = float(acc_summ_df[acc_summ_df['Tag'] == 'UnrealizedPnL']['Value'].values[0]) # USD



