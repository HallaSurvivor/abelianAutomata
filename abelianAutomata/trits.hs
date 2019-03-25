module Trits where

enum = enum' 0
  where
    enum' _ [] = []
    enum' n (x:xs) = (x,n):(enum' (n+1) xs)

trits 0 = [[]]
trits n = map (1:) t ++ map (0:) t ++ map (-1:) t
  where
    t = trits (n-1)

eval ts = sum $ map (\(t,n) -> t * 2^n) $ enum ts
