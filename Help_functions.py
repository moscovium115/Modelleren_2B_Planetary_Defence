import numpy as np

def calculate_richardson(max_arr):
    print("lengtes", len(max_arr[0]))
    print("lengtes", len(max_arr[1][0::2]))
    print("lengtes", len(max_arr[2][0::4]))

    # max_arr[1]=max_arr[1][0::2]
    # max_arr[2]=max_arr[2][0::4]
    # rich_arr=(max_arr[1]-max_arr[0])/ (max_arr[2]-max_arr[1])

    #
    # # we houden alleen de positieve waarden
    # rich_arr=rich_arr[rich_arr>0]
    # print("richardson error:",rich_arr,"fout orde:", np.log(np.nanmean(rich_arr))/np.log(2))

    arr_4H=max_arr[0]
    arr_2H=max_arr[1][0::2]
    arr_H=max_arr[2][0::4]
    print(" x debuggggg ", max_arr[0])

    #if all array lengths are equal, we can directly calculate the richardson error
    if len(arr_4H)==len(arr_2H):
        if len(arr_2H)==len(arr_H):
            rich_arr=(arr_2H-arr_4H)/ (arr_H-arr_2H)
            # we houden alleen de positieve waarden
            rich_arr=rich_arr[rich_arr>0]
            print("richardson error:",rich_arr,"fout orde:", np.log(np.nanmean(rich_arr))/np.log(2))
        else:
            return arr_4H, arr_2H

    #else we need to remove some values to make them equal in length
    else:
        if len(arr_4H)>len(arr_2H):
            arr_4H=np.delete(arr_4H, -1)
        else:
            arr_2H=np.delete(arr_2H, -1)


    # if (len(max_arr[1][1::2])>len(max_arr[0])):
    #     max_arr[1]=np.delete(max_arr[1][1::2],-1)
    # else:
    #     max_arr[1]=max_arr[1][1::2]
    # if (len(max_arr[2][3::4])>len(max_arr[0])):
    #     max_arr[2]=np.delete(max_arr[2][3::4],-1)
    # else:
    #     max_arr[2]=max_arr[2][3::4]
    #
    #
    # rich_arr=(max_arr[1]-max_arr[0])/ (max_arr[2]-max_arr[1])
    # print("richardson:",np.nanmean(rich_arr), "fout orde:", np.log(np.nanmean(rich_arr))/np.log(2))


