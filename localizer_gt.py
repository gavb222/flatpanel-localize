import torch
import math

#TODO: Ensure odd output
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
#x, y = center of the gaussian
#spread = spread of the gaussian
#NB that x-spread > 0, y-spread > 0, x+spread < panel_x, y+spread < panel_y
def produce_freq_response(panel_x, panel_y, n_freq, x, y, spread):

    #TODO: change these to return errors
    if x-spread < 0:
        return -1
    elif y-spread < 0:
        return -1
    elif x+spread > panel_x:
        return -1
    elif y+spread > panel_y:
        return -1

    response = gaussian(spread)

    #response.size is (dim,dim)
    n_gaussian_elems = response.size()[0]

    #pad response with zeros until its the size we want
    n = math.floor(n_gaussian_elems/2)
    #pad of x starting from 0,
    pad_left = torch.zeros((x-n),n_gaussian_elems)
    pad_right = torch.zeros((panel_x-(x+n))-1,n_gaussian_elems)

    pad_top = torch.zeros(panel_x,(y-n))
    pad_bottom = torch.zeros(panel_x,(panel_y-(y+n))-1)

    response = torch.cat((pad_left,response), dim=0)
    response = torch.cat((response,pad_right), dim=0)
    response = torch.cat((pad_top,response), dim=1)
    response = torch.cat((response,pad_bottom), dim=1)

    response = response.expand(n_freq,panel_x,panel_y)
    return response

gt = produce_freq_response(10,10,2,9,9,5)
print(gt)
