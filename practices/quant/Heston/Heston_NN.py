import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import grad


class Net(nn.Module):

    def __init__(self,num_input=7, num_neurons=128, device=None):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(num_input, num_neurons) 
        self.fc2 = nn.Linear(num_neurons, num_neurons)
#        self.fc3 = nn.Linear(num_neurons, num_neurons)
#        self.fc4 = nn.Linear(num_neurons, num_neurons)
#        self.fc5 = nn.Linear(num_neurons, num_neurons)
#        self.fc6 = nn.Linear(num_neurons, num_neurons)
        self.fc7 = nn.Linear(num_neurons, 1)

        if device is None:
            self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device
        self.to(self.device)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = torch.tanh(self.fc2(x))
#        x = torch.tanh(self.fc3(x))
#        x = F.relu(self.fc3(x))
#        x = F.relu(self.fc4(x))
#        x = F.relu(self.fc5(x))
#        x = F.relu(self.fc6(x))
        x = F.relu(self.fc7(x))
#        x = torch.sigmoid(self.fc7(x))
        return x
    
    def get_loss(self, x, y_train, T_loc, k_loc, criterion, arbitrage_weight, l2_weight, show_log=False):
        y_hat = self.forward(x)
        loss = criterion(y_hat, y_train)
        l2_reg = torch.tensor(0.).to(self.device)

        # Penalization (loop over each data)
        ## calendar arbitrage
        calendar_arbi_count = 0
        calendar_loss = torch.tensor(0.).to(self.device)
        butterfly_loss = torch.tensor(0.).to(self.device)
        for elem in x:
            y= self.forward(elem)
            dydx = grad(y, elem, create_graph = True)[0]
            dydT = dydx[T_loc]
            if dydT < 0.0:
                calendar_arbi_count += 1
                calendar_loss += (torch.exp(-dydT) * arbitrage_weight[0])
        loss += calendar_loss
        ## butterfly arbitrage
        butterfly_count = 0
        dydk = dydx[k_loc] # dCdk w.r.t. log-strike 
        d2ydk2 = grad(dydx[k_loc],elem, create_graph = True)[0][k_loc] # d2Cdk2 w.r.t. log-strike
        # butterfly arbitrage: d2CdK2 = 1/K^2(d2Cdk2 - dCdk) > 0, d2Cdk2 > dCdk
        # dCdK = 1/K*dcdk
        # d2CdK2 = 1/e^2k * (d2Cdk2 - dCdk)
        butter_ineq = (d2ydk2 - dydk)/torch.exp(2.0*elem[k_loc])
        if butter_ineq < 0.0: # violation of butterfly arbitrage
            butterfly_count += 1
            butterfly_loss += torch.exp(-butter_ineq) * arbitrage_weight[1]
        loss += butterfly_loss

        # regularizations
        for param in self.parameters():
            l2_reg += torch.norm(param)
            loss += l2_weight * l2_reg
            
        if show_log:
            if calendar_arbi_count > 0: print('calendar ',calendar_arbi_count,' ',calendar_loss)
            if butterfly_count > 0: print('butterfly ',butterfly_count,' ',butterfly_loss)
            print(f'Loss: {loss}')
            
        return loss
        
    def get_loss_test(self, x, y_train, T_loc, k_loc, criterion, arbitrage_weight, l2_weight, show_log=False):
        y_hat = self.forward(x)
        loss = criterion(y_hat, y_train)
        l2_reg = torch.tensor(0.).to(self.device)

        calendar_arbi_count = 0
        butterfly_count = 0
        # Penalization (loop over each data)
        if sum(arbitrage_weight) > 0.0:            
            calendar_loss = torch.tensor(0.).to(self.device)
            butterfly_loss = torch.tensor(0.).to(self.device)
            for elem in x:
                ## calendar arbitrage
                y= self.forward(elem)
                dydx = grad(y, elem, create_graph = True)[0]
                if arbitrage_weight[0] > 0.0:
                    dydT = dydx[T_loc]
                    if dydT < 0.0:
                        calendar_arbi_count += 1
                        calendar_loss += torch.exp(-dydT)
                if arbitrage_weight[1] > 0.0:
                    ## butterfly arbitrage
                    dydk = dydx[k_loc] # dCdk w.r.t. log-strike 
                    d2ydk2 = grad(dydk,elem, create_graph = True)[0][k_loc] # d2Cdk2 w.r.t. log-strike
                    # butterfly arbitrage condition: d2CdK2 = 1/K^2(d2Cdk2 - dCdk) > 0, d2Cdk2 > dCdk
                    # dCdK = 1/K*dcdk
                    # d2CdK2 = 1/e^2k * (d2Cdk2 - dCdk)
                    butter_ineq = torch.sub(d2ydk2, dydk)
                    if butter_ineq < 0.0: # violation of butterfly arbitrage
                        butter_ineq = torch.div(butter_ineq, torch.exp(2.0*elem[k_loc]))
                        butterfly_count += 1
                        butterfly_loss += torch.exp(-butter_ineq)
            calendar_loss *= arbitrage_weight[0]
            butterfly_loss *= arbitrage_weight[1]
            loss += calendar_loss
            loss += butterfly_loss

        # regularizations
        if l2_reg > 0.0:            
            for param in self.parameters():
                torch.add(l2_reg, torch.norm(param))
            loss += torch.mul(l2_reg, l2_weight)
            
        if show_log:
            if calendar_arbi_count > 0: print('calendar ',calendar_arbi_count,' ',calendar_loss)
            if butterfly_count > 0: print('butterfly ',butterfly_count,' ',butterfly_loss)
            print(f'Loss: {loss}')
            
        return loss

    def get_loss_fd(self, x, y_train, T_loc, k_loc, criterion, arbitrage_weight, l2_weight, show_log=False):
        y_hat = self.forward(x)
        loss = criterion(y_hat, y_train)
        l2_reg = torch.tensor(0.).to(self.device)
        incr = 0.05 # stencil ratio

        # calendar arbitrage
        xp = x.clone().detach()        
        xp[:,T_loc] += incr
        yp = self.forward(xp)
        dydT = (yp -y_hat) / incr
        mask = dydT < 0.0
        if torch.sum(mask) > 0.0:
            calendar_loss = torch.sum(torch.exp(dydT[dydT<0.0]))
            loss += calendar_loss*arbitrage_weight[0]

        # butterfly arbitrage
        x_tmp = x.clone().detach()        
        x_tmp[:,k_loc] += incr
        yp = self.forward(x_tmp)
        x_tmp[:,k_loc] -= 2.0*incr
        ym = self.forward(x_tmp)
        dydk = (yp - ym) / (2.0*incr)
        d2ydk2 = (yp + ym - 2.0*y_hat) / (incr**2)
        butter_ineq = torch.sub(d2ydk2, dydk)
        butter_ineq = butter_ineq.reshape(-1)
        mask = butter_ineq < 0.0
#        print(mask)
        if torch.sum(mask) > 0.0:
#            print(butter_ineq)
#            print(x)
#            print(butter_ineq[mask])
#            return mask, x
#            print(x[mask])
            butterfly_loss = torch.sum(torch.exp(-butter_ineq[mask] / torch.exp(2.0*x[mask,k_loc])))        
            loss += butterfly_loss*arbitrage_weight[1]


        # regularizations
        if l2_reg > 0.0:            
            for param in self.parameters():
                torch.add(l2_reg, torch.norm(param))
            loss += torch.mul(l2_reg, l2_weight)

#        if show_log:
#            if calendar_arbi_count > 0: print('calendar ',calendar_arbi_count,' ',calendar_loss)
#            if butterfly_count > 0: print('butterfly ',butterfly_count,' ',butterfly_loss)
#            print(f'Loss: {loss}')

        return loss

    
def weights_init(m):
    if isinstance(m, nn.Linear):
        torch.nn.init.xavier_uniform_(m.weight)
#        torch.nn.init.normal_(m.weight)
#        xavier(m.weight.data)
#        xavier(m.bias.data)
