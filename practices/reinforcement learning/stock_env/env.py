import random
import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import pandas as pd

class TradingSPYEnv(gym.Env):
    """
    SPY (S&P500) trading environment.
  
    State: [[short, neutral, long], portfolio value]
      - The states are 
  
  
    Action: sell (0), hold (1), and buy (2)
      - I prescribe a very simple policy
      - when selling, sell all the shares
      - when buying, buy as many as cash in hand allows
    """
    def __init__(self, train_data_path='historySPY.csv', sma_len=[5], init_invest=10000):
        train_data = pd.read_csv(train_data_path, index_col = False, parse_dates= ['Date'])
        self.stock_price_history = train_data 
        self.current_step = 0 # The step in the data
        self.iteration = 0 # the iteration step in an episode
        self.init_invest = init_invest
        self.accumulated_profit = 0.0

        feature_dict = {'Date': self.stock_price_history['Date'],
                    'State': np.zeros(self.stock_price_history.shape[0], dtype=int),
                    'accumulated_profit': np.zeros(self.stock_price_history.shape[0], dtype=float) ,
                    'portfolio_value': np.zeros(self.stock_price_history.shape[0], dtype=float),
                    'Close': self.stock_price_history['Close']
                    }
    
        # feature engineering. Put values like sma
        feature = 'Close'
        if isinstance([],list):
            for sma in sma_len:
                feature_dict[feature+'_'+str(sma)] = self.stock_price_history[feature].rolling(sma).mean()
#            self.stock_price_history[feature+'_'+str(sma)] = self.stock_price_history[feature].rolling(sma).mean()
                    
        self.stock_price_history.dropna(axis=0,inplace=True)
        self.stock_price_history.reset_index(drop=True,inplace=True)
        self.features = pd.DataFrame(feature_dict).dropna(axis=0)
        self.features.reset_index(drop=True,inplace=True)
        self.features.portfolio_value.loc[self.current_step] = self.init_invest
    
        self.end_step = self.features.shape[0]
    
        # action space
        # 0: short, 1: neutral, 2: long
        self.action_space = spaces.Discrete(3)
    
        # observation space
        # This contains features to make decisions
        self.observation_space = spaces.Box(low=0, high=np.inf, shape=(self.features.columns.shape[0] -1,), dtype=np.float16)
    
    
    def _get_observation(self):
        # return features at current step
        tmp = self.features.drop(columns=['Date']).loc[self.current_step].to_numpy()
        # state, portfolio_value, Close, smas
        return tmp

    def reset(self):
        self.iteration = 0 
        self.features['State'] = 1 # State:1 means market neutral
        self.features['portfolio_value'] = 0.0       
        self.features['accumulated_profit'] = 0.0
        
        # Set the current step to a random point within the data frame
        self.current_step = random.randint(
            0, int(self.features.shape[0] * 0.9))
        self.features['portfolio_value'].loc[self.current_step] = self.init_invest

        return self._get_observation()

    """
    Compute what happens next step
    """
    def step(self, action):
        next_step = self.current_step + 1
        prev_step = self.current_step - 1
        if next_step == self.end_step: 
            # At the end, we have nothing to do
            done = True
            return None, None, done, None
        
        col_name = 'Close'
        features = self.features
        portfolio_value = self.features.portfolio_value        

        done = False
        # reward 
        r_t = 0.0
        if (self.iteration > 1) and (self.current_step is not self.end_step): # Exclude the very first step
            # difference in portfolio value 
            r_t = portfolio_value.loc[self.current_step] - portfolio_value.loc[prev_step]
            

        features['accumulated_profit'].loc[self.current_step] = features['accumulated_profit'].loc[prev_step] + r_t    
            
        # Current state is set
        self.features.State.loc[self.current_step] = action
    
        # Compute next step
        # compute portfolio value at next step
        if action == 0: # shorting
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step] * features[col_name].loc[self.current_step] / features[col_name].loc[next_step]
        elif action == 1: # market-neutral position (100% cash)  
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step]
        elif action == 2: # longing
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step] * features[col_name].loc[next_step] / features[col_name].loc[self.current_step]
        else:
            raise TypeError("Action is out of the space")
        self.features.State.loc[next_step] = action
    
        self.current_step += 1 
        self.iteration += 1
        s_prime = self._get_observation() # state at t+1
    
        return s_prime, r_t, done, None
