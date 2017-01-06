module Dna.Search (
    dnaFromSepll, composingKMer,
    initDeBrujinGraphFromKMers,
    initDeBrujinGraphFromDna, fibrabbit) where

import Dna.Graph
import Dna.Nucleotides
import Data.List
import Control.Applicative
import Data.Maybe
import Control.Monad
import Numeric.LinearAlgebra.Devel

composingKMer ::  Int -> Dna -> DnaList
composingKMer n dna
    | lengthDna dna >= n = takeDna n dna : composingKMer n (tailDna dna)
    | otherwise = []

same one two = one `isSubsequenceOf` two && two `isSubsequenceOf` one

mergeHeadTail [] r = r
mergeHeadTail (lx:lxs) r@(rx:rxs)
    | lx == rx = r
    | otherwise = lx: mergeHeadTail lxs r

isConnected [] _ = False
isConnected _ [] = False
isConnected (_:lxs) r = lxs `isPrefixOf` r


dnaFromSepll :: DnaList -> Dna
dnaFromSepll [] = nullDna
dnaFromSepll [dna] = dna
dnaFromSepll (x:xs) = dnaSpellStep 1 x xs
    where
        dnaSpellStep dropStep growDna dnas = case find (isDnaPrefixOf $ dropDna dropStep growDna) dnas of
            Just newDna -> dnaSpellStep (dropStep+1) (takeDna dropStep growDna +:+ newDna) dnas
            Nothing -> growDna


initDeBrujinGraphFromDna :: Int -> Dna -> Graph Dna
initDeBrujinGraphFromDna n dnas = initOverlaGraph nodes $ \ m findIndex ->
        forM_ adjacentNodes $ \ (from, to) ->
            modifyMatrix m (findIndex from) (findIndex to) (1+)
    where
        nodes = composingKMer n dnas
        adjacentNodes = zip (init nodes) (tail nodes)

initDeBrujinGraphFromKMers :: DnaList -> Graph Dna
initDeBrujinGraphFromKMers kmers = initOverlaGraph nodes $ \m findIdxNode ->
    forM_ nodes $ \from ->
        forM_ (map tailDna (filter (from `isDnaPrefixOf`) kmers)) $ \to ->
            modifyMatrix m (findIdxNode from) (findIdxNode to) (1+)
        where nodes = nub $ map (liftDna init) kmers ++ map tailDna kmers


fibrabbit 1 _ = 1
fibrabbit 2 _ = 1
fibrabbit n k = fibrabbit (n-1) k + (k * fibrabbit (n-2) k)
