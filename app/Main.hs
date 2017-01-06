module Main (main) where

import Dna.Search
import Data.Maybe
import Dna.Nucleotides
import Dna.MetaInfo
import Control.Monad
import Dna.Graph
import Data.Maybe
import Debug.Hoed.Pure

dnaSpell = do
    dnas <- read <$> readFile ""
    print $ dnaFromSepll dnas

readDnasFromFile :: String -> IO DnaList
readDnasFromFile = readDnas . readFile

readDna :: IO String -> IO Dna
readDna io = read <$> io

readDnas :: IO String -> IO DnaList
readDnas io = read <$> io

readInt :: IO String -> IO Int
readInt io = read <$> io

overlapMain = do
    dnas <- readDnasFromFile "rsamples/osalind_ba3c.txt"
    let out  = show $ overlapGraphFromNode dnas
    writeFile "samples/overlap_graphout.txt" out
    print "Finished"

debrujin1 = do
    let ioFile = lines <$> readFile "rsamples/osalind_ba3d.txt"
    len <- readInt $ liftM2 (!!) ioFile (return 0)
    dna <- readDna $ liftM2 (!!) ioFile (return 1)
    let out = show $  initDeBrujinGraphFromDna (len-1) dna
    writeFile "samples/debrujin_graphout.txt" out
    print "finished"

debrujin2 = do
    dnas <- readDnasFromFile "samples/debrunjin_list_sample.txt"
    let out = show $ initDeBrujinGraphFromKMers dnas
    writeFile "samples/debrujin_list_graphout.txt" out
    print "Finished"

complement = do
    dna <- readLn
    print $ dnaComplement dna

canBeBalanceTest :: IO ()
canBeBalanceTest = do
    inGraph <- readFile "samples/graph1_ba3g.txt"
    let graph = (graphFromConnection $ readAllConnection inGraph) :: Graph Int
    print $ canBeBalance graph

balanceTest :: IO ()
balanceTest = do
    inGraph <- readFile "samples/graph1_ba3g.txt"
    let graph = (graphFromConnection $ readAllConnection inGraph) :: Graph Int
    print $ balance graph

euclidianCycleTest = do
    inGraph <- readFile "samples/graph1_ba3g.txt"
    let graph = (graphFromConnection $ readAllConnection inGraph) :: Graph Int
    let cycle = euclidianCycle graph
    let out = liftM2 (++) (fmap (show.fst.head) cycle) (fmap (concatMap (("->" ++) . show . snd)) cycle)
    writeFile "samples/sample_graphout.txt" (show (fromMaybe "" out))

main = do
    fastaList <- parseFastaFile "samples/fasta_ex.fas"
    return $ fmap highestGC fastaList
