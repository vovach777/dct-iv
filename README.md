# port from python

```python
import numpy as np

def dctIV(audio_signal):
    window_length = len(audio_signal)
    audio_dct = np.zeros((8*window_length))

    audio_dct[1:2*window_length:2] = audio_signal
    audio_dct[2*window_length+1:4*window_length:2] = -np.flip(audio_signal)
    audio_dct[4*window_length+1:6*window_length:2] = -audio_signal
    audio_dct[6*window_length+1:8*window_length:2] = np.flip(audio_signal)
    #print(repr(audio_dct))
    audio_dct = np.fft.fft(audio_dct, axis=-1, norm=None)
    #print(repr(audio_dct))
    audio_dct = np.real(audio_dct[1:2*window_length:2]) / 4

    audio_dct *= np.sqrt(2/window_length)

    return audio_dct


np.set_printoptions(linewidth=100,threshold=99999999)
# Пример использования
audio_signal = np.array([ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])  # Пример входных данных
result = dctIV(audio_signal)
print(result)
result = dctIV(result)
print(result)
```
