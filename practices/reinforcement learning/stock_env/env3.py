import random
import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import pandas as pd

class TradingSPYEnv(gym.Env):
    """
    SPY (S&P500) trading environment.
  
    The features are boolean. Signal of breakout smas

  
    Action: sell (0), hold (1), and buy (2)
      - I prescribe a very simple policy
      - when selling, sell all the shares
      - when buying, buy as many as cash in hand allows
    """
    def __init__(self, train_data_path='historySPY.csv', sma_len=[5], init_invest=10000, learning_rate=0.0002, gamma=0.98,
                normalize_price = True, mode = 'train', train_test_split = 0.9,len_history_state=0):
        train_data = pd.read_csv(train_data_path, index_col = False, parse_dates= ['Date'])
        self.stock_price_history = train_data 
        self.max_sma_len = max(sma_len)
        self.len_history_state = len_history_state
        self.index_offset = self.max_sma_len + len_history_state
        self.current_step = self.max_sma_len+self.index_offset # The step in the data
        self.iteration = 0 # the iteration step in an episode
        self.init_invest = init_invest
        self.accumulated_profit = 0.0
        self.normalize_price = normalize_price

        feature_dict = {'Date': self.stock_price_history['Date'],
                    'State': np.zeros(self.stock_price_history.shape[0], dtype=int),
                    'accumulated_profit': np.zeros(self.stock_price_history.shape[0], dtype=float), 
                    'portfolio_value': np.zeros(self.stock_price_history.shape[0], dtype=float),
                    'Close': self.stock_price_history['Close']
                    }
    
        # feature engineering. Put values like sma
#        if sma_len([],list):
            
#            for sma in sma_len:
#                feature_dict[feature+'_'+str(sma)] = self.stock_price_history[feature].rolling(sma).mean()
#            self.stock_price_history[feature+'_'+str(sma)] = self.stock_price_history[feature].rolling(sma).mean()
                    
        self.stock_price_history.dropna(axis=0,inplace=True)
        self.stock_price_history.reset_index(drop=True,inplace=True)

        self.features = pd.DataFrame(feature_dict)
        if isinstance(sma_len,list):
            self._set_sma(sma_len)
            self._set_breakout_sma(sma_len)
            self._set_daily_return()
    
        train_test_split_index = int(self.features.shape[0] * train_test_split)
        if mode == 'train':
            self.end_step = train_test_split_index
        elif mode == 'test':
            self.features.shape[0]
            self.current_step = train_test_split_index
            self.end_step = self.features.shape[0]

        # Set up data and features
        self.reset(current_step = self.current_step)
                    
        # action space
        # 0: short, 1: neutral, 2: long
        self.action_space = spaces.Discrete(3)
    
        # observation space
        # This contains features to make decisions
        self.observation_space = spaces.Box(low=0, high=np.inf, shape=(self.features.columns.shape[0] -1,), dtype=np.float16)
    
    def _set_sma(self, sma_len):
        feature = 'Close'
        for sma in sma_len:
            sma_name = feature+'_'+str(sma)
            self.stock_price_history[sma_name] = self.stock_price_history[feature].rolling(sma).mean()
            self.features[sma_name] = self.stock_price_history[feature].rolling(sma).mean()
        self.features = self.features.dropna(axis=0)        
            
    def _set_breakout_sma(self,sma_len):
        feature = 'Close'
        features = self.features
        for sma in sma_len:
            sma_name = feature+'_'+str(sma)
            breakout_name = 'breakout_' + str(sma)
            features[breakout_name] = features[feature] >= features[sma_name]
            features.drop(columns=[sma_name],inplace=True)
        
    def _set_daily_return(self):
        self.features['daily_return'] = self.stock_price_history['Close'].pct_change()
        
    def _get_observation(self):
        # return features at current step
#        print(self.features)
#        print(self.current_step)
#        observation = self.features.drop(columns=['Date']).loc[self.current_step].to_numpy('float32')

        # Only return things I want to
        observation = []
        observation.append(self.features['State'].loc[self.current_step-self.len_history_state+1:self.current_step].to_numpy('float32'))
        observation.append(self.features['accumulated_profit'].loc[self.current_step-self.len_history_state+1:self.current_step].to_numpy('float32'))
        observation.append(self.features['daily_return'].loc[self.current_step-self.len_history_state+1:self.current_step].to_numpy('float32'))
        # breakout signals
        for col in self.features.columns:
            if 'breakout_' in col:
                observation.append(self.features[col].loc[self.current_step-self.len_history_state+1:self.current_step].values)
            
        return observation

    def reset(self, current_step = None):
        self.iteration = 0 
        self.features['State'] = 1 # State:1 means market neutral
        self.features['portfolio_value'] = 0.0       
        self.features['accumulated_profit'] = 0.0
        
        # Set the current step to a random point within the data frame
        if current_step is not None:
            self.current_step = current_step
        else:
            self.current_step = random.randint(self.index_offset, int(self.features.shape[0] * 0.9))
            
        self.features['portfolio_value'].loc[self.current_step] = self.init_invest
        
        # normalize Close price
        price = self.stock_price_history['Close'].loc[self.current_step]
        self.features['Close'].loc[self.current_step:self.end_step] = self.stock_price_history['Close'].loc[self.current_step:self.end_step] / price
        
#        if self.normalize_price:
#            price = self.stock_price_history['Close'].loc[self.current_step]
#            self.features[].loc[self.current_step:self.end_step] = self.stock_price_history[col].loc[self.current_step:self.end_step] / price
#            for col in self.features.columns:
#                if 'Close' in col:
#                    self.features[col].loc[self.current_step:self.end_step] = self.stock_price_history[col].loc[self.current_step:self.end_step] / price

        return self._get_observation()

    """
    Compute what happens next step
    """
    def step(self, action):
        next_step = self.current_step + 1
        prev_step = self.current_step - 1
        
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

        if next_step == self.end_step: 
            # At the end, we have nothing to do
            done = True
            self.accumulated_profit = features['accumulated_profit'].loc[self.current_step]
            long_return = self.features['Close'].loc[self.current_step] / self.features['Close'].loc[self.current_step-self.iteration]
            return None, None, done, {'profit_iteration': self.accumulated_profit/self.iteration, 'iterations': self.iteration, 
                                      'long_return': long_return}
        
        
        # Current state is set
        self.features.State.loc[self.current_step] = action
    
        # Compute next step
        # compute portfolio value at next step
        stock_price_history = self.stock_price_history
        if action == 0: # shorting
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step] * stock_price_history[col_name].loc[self.current_step] / stock_price_history[col_name].loc[next_step]
        elif action == 1: # market-neutral position (100% cash)  
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step]
        elif action == 2: # longing
            portfolio_value.loc[next_step] = portfolio_value.loc[self.current_step] * stock_price_history[col_name].loc[next_step] / stock_price_history[col_name].loc[self.current_step]
        else:
            raise TypeError("Action is out of the space")
        self.features.State.loc[next_step] = action
    
        self.current_step += 1 
        self.iteration += 1
        s_prime = self._get_observation() # state at t+1
    
        return s_prime, r_t, done, None
