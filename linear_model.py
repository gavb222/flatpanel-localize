import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import time
import random
import matlab.engine

def gaussian(spread):
    #spread controls the size of the array
    linspace = torch.linspace(-2.5,2.5,spread)
    # gaussian = e^((-x)^2/2) when standard dev is 1 and height is 1
    linspace = torch.exp(-1 * torch.div(torch.pow(linspace,2),2))

    out_x = linspace.expand(spread,spread)
    out_y = out_x.permute(1,0)
    out_gaussian = out_x * out_y
    return out_gaussian

#panel_x, panel_y = panel dimensions
#n_freq = n frequency bins
#x, y = top left of gaussian
#spread = spread of the gaussian
#NB that x-spread > 0, y-spread > 0, x+spread < panel_x, y+spread < panel_y
def produce_freq_response(panel_x, panel_y, n_freq, x, y, spread, expand_dim=False):

    #TODO: change these to return errors
    if x+spread > panel_x-1:
        return torch.ones(panel_x,panel_y)*-1
    elif y+spread > panel_y-1:
        return torch.ones(panel_x,panel_y)*-1

    response = gaussian(spread)

    #response.size is (dim,dim)
    #n_gaussian_elems = response.size()[0]

    #pad response with zeros until its the size we want
    #n = math.floor(n_gaussian_elems/2)

    #pad of x starting from 0,
    #pad_left = torch.zeros((x-n),n_gaussian_elems)
    #pad_right = torch.zeros((panel_x-(x+n))-1,n_gaussian_elems)

    #pad_top = torch.zeros(panel_x,(y-n))
    #pad_bottom = torch.zeros(panel_x,(panel_y-(y+n))-1)

    #response = torch.cat((pad_left,response), dim=0)
    #response = torch.cat((response,pad_right), dim=0)
    #response = torch.cat((pad_top,response), dim=1)
    #response = torch.cat((response,pad_bottom), dim=1)

    out_array = torch.zeros(panel_x,panel_y)
    out_array[x:x+spread,y:y+spread] = response

    if expand_dim:
        out_array = out_array.expand(n_freq,panel_x,panel_y)

    return out_array

class Conv_Block(nn.Module):
    def __init__(self, input_size, output_size, kernel_size=4, stride=2, padding=1,activation=True):
        super(Conv_Block, self).__init__()
        self.conv = nn.Conv2d(input_size, output_size, kernel_size, stride, padding)
        self.activation = activation

    def forward(self, x):
        if self.activation:
            out = self.conv(F.relu(x))
        else:
            out = self.conv(x)

        return out

class Conv_Net(nn.Module):
    def __init__(self, input_channels, n_filters, output_channels):
        super(Conv_Net, self).__init__()

        self.conv1 = Conv_Block(input_channels, n_filters, activation=False)
        self.conv2 = Conv_Block(n_filters, n_filters * 2)
        self.conv3 = Conv_Block(n_filters * 2, n_filters * 4)
        self.conv4 = Conv_Block(n_filters * 4, n_filters * 8, stride=1)
        self.conv5 = Conv_Block(n_filters * 8, output_channels, stride=1)
        self.classifier = nn.Linear(384,24)

    def forward(self,x):
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)
        x = self.conv5(x)
        x = x.view(-1)
        x = self.classifier(x)
        #what size is this?
        out = torch.nn.Sigmoid()(x)
        out = torch.reshape(out,(4,6))
        return out

model = Conv_Net(1,16,24)
model.cuda()
model.train()

loss_fn = nn.MSELoss
criterion = torch.optim.Adam(model.parameters(), lr = .0001, betas = (.5,.999))

keep_training = True
epoch_counter = 0

panel_x = 50
panel_y = 50

eng = matlab.engine.start_matlab()

#make a panel
driver_locations = torch.tensor((0.25, 0.25, 0.75, 0.75, 0.25, 0.75, 0.75, 0.25)).view(4,2)
Lx = 0.3
Ly = 0.5


while keep_training:
    epoch_counter = epoch_counter + 1
    time_start = time.time()
    gt = torch.ones(panel_x,panel_y)*-1

    model.zero_grad()

    #random init starting conditions
    while gt[0,0] == -1:
        #returns -1 for invalid configuration
        gt = produce_freq_response(panel_x,panel_y,1,random.randint(1,panel_x-1),random.randint(1,panel_y-1),random.randint(3,15))

    coefs = model(gt.unsqueeze(0).unsqueeze(0).cuda())

    print(coefs.size())
    #very possible that the interpreter doesnt like torch tensors, might have to go numpy with this
    response1, frequencies = eng.get_biquad_response(coefs[0,:].cpu().detach().numpy(),44100,nargout = 2)
    response2, temp = eng.get_biquad_response(coefs[1,:].cpu().detach().numpy(),44100,nargout = 2)
    response3, temp = eng.get_biquad_response(coefs[2,:].cpu().detach().numpy(),44100,nargout = 2)
    response4, temp = eng.get_biquad_response(coefs[3,:].cpu().detach().numpy(),44100,nargout = 2)

    responses = torch.stack((response1,response2,response3,response4),dim=-1)

    matlab_panel = eng.Clamped_Panel[driver_locations,responses,frequencies,Lx,Ly]

    matlab_out = eng.matlab_panel.view_total_scan(200,0)

    loss = loss_fn(matlab_out,gt)
    criterion.step()

    print("holy moly!")
