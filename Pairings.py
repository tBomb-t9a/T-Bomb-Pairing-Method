import numpy as np
import pandas as pd
from itertools import combinations
from scipy.optimize import linear_sum_assignment


def readMatrix(ours, theirs, matrix):
    return(ours, theirs, pd.DataFrame(matrix, index = ours, columns = theirs))

def predictedOutcome(ourPlayer, theirPlayer, matrix):
    prediction = matrix[theirPlayer][ourPlayer]
    return prediction

def pairingCombos(players):
    return (list(combinations(players, 2)))

def round2Pick(ourFirstPick, theirFirstPick, matrix):
    ourPlayers = list(matrix.index.values)
    theirPlayers = list(matrix.columns.values)

    ourPlayers.remove(ourFirstPick)
    theirPlayers.remove(theirFirstPick)

    return(pairingCombos(ourPlayers), pairingCombos(theirPlayers))

def notPicked(ourFirstPick, theirFirstPick, ourRound2Picks, theirRound2Picks, matrix):
    ourPlayers = list(matrix.index.values)
    theirPlayers = list(matrix.columns.values)

    ourPlayers.remove(ourFirstPick)
    theirPlayers.remove(theirFirstPick)

    for player in ourRound2Picks:
        ourPlayers.remove(player)

    for player in theirRound2Picks:
        theirPlayers.remove(player)

    return((ourPlayers[0], theirPlayers[0]))

def finalScores(ourFirstPick, theirFirstPick, ourRound2Picks, theirRound2Picks, matrix):
    playersNotPicked = notPicked(ourFirstPick, theirFirstPick, ourRound2Picks, theirRound2Picks, matrix)

    ourChoice = (playersNotPicked[0], ourFirstPick)
    theirChoice = (playersNotPicked[1], theirFirstPick)

    ourBestPairing = pairScoresUs(ourChoice, theirRound2Picks, matrix)
    theirBestPairing = pairScoresThem(ourRound2Picks, theirChoice, matrix)

    return ([ourBestPairing[0], ourBestPairing[1], theirBestPairing[0], theirBestPairing[1],  ourBestPairing[2] + theirBestPairing[2]])

def pairScoresUs(ourTwoPlayers, theirTwoPlayers, matrix):
    score1 = predictedOutcome(ourTwoPlayers[0], theirTwoPlayers[0], matrix) + predictedOutcome(ourTwoPlayers[1], theirTwoPlayers[1], matrix)
    score2 = predictedOutcome(ourTwoPlayers[1], theirTwoPlayers[0], matrix) + predictedOutcome(ourTwoPlayers[0], theirTwoPlayers[1], matrix)
    if score1 >= score2:
        return([(ourTwoPlayers[0], theirTwoPlayers[0]), (ourTwoPlayers[1], theirTwoPlayers[1]), score1])
    else:
        return([(ourTwoPlayers[1], theirTwoPlayers[0]), (ourTwoPlayers[0], theirTwoPlayers[1]), score2])

def pairScoresThem(ourTwoPlayers, theirTwoPlayers, matrix):
    score1 = predictedOutcome(ourTwoPlayers[0], theirTwoPlayers[0], matrix) + predictedOutcome(ourTwoPlayers[1], theirTwoPlayers[1], matrix)
    score2 = predictedOutcome(ourTwoPlayers[1], theirTwoPlayers[0], matrix) + predictedOutcome(ourTwoPlayers[0], theirTwoPlayers[1], matrix)
    if score1 <= score2:
        return([(ourTwoPlayers[0], theirTwoPlayers[0]), (ourTwoPlayers[1], theirTwoPlayers[1]), score1])
    else:
        return([(ourTwoPlayers[1], theirTwoPlayers[0]), (ourTwoPlayers[0], theirTwoPlayers[1]), score2])

def averageScore(ourPlayer, theirPlayer, matrix):
    possiblePairings = round2Pick(ourPlayer, theirPlayer, matrix)
    results = []

    for ourPair in possiblePairings[0]:
        score = 0
        for theirPair in possiblePairings[1]:
            score += finalScores(ourPlayer, theirPlayer, ourPair, theirPair, matrix)[-1]
        results.append((ourPair, score/len(possiblePairings[1])))

    return(results)

def firstPickAvg(ourPlayer, theirPlayer, matrix):
    res = averageScore(ourPlayer, theirPlayer, matrix)
    return (res[0][1] + res[1][1] + res[2][1])/3 # Change this, it sucks

def initialScores(ours, theirs, matrix):
    for us in ours:
        score = 0
        for them in theirs:
            score += firstPickAvg(us, them, matrix)
        print("If you pick {} first, the average score is {}".format(us, score/len(theirs)))
    return

def getFirstPlayer(choices):
    choice = ""
    while choice not in choices:
        choice = input("Choose one of [%s]:" % ", ".join(choices))
    return(choice)

def printAverageScores(ourFirstPick, theirFirstPick, matrix):
    avgsc = averageScore(ourFirstPick, theirFirstPick, matrix)
    print("For round 2 the average points our picks will get is: ")
    for score in avgsc:
        print("{} and {}: {} points".format(score[0][0], score[0][1], score[1]))
    return

def main():
    ours = ["BH", "DL", "DE", "DH"]
    theirs = ["EoS", "HE", "ID", "KoE"]
    scoreArray = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]).reshape((4,4))
    theirScores = np.ones(16).reshape((4,4))*6 - scoreArray
    df = pd.DataFrame(scoreArray, index = ours, columns = theirs)
    # ours, theirs, df = readMatrix(ours, theirs, A)
    print(df)
    row_ind, col_ind = linear_sum_assignment(theirScores)
    matchingAlgRes = scoreArray[row_ind, col_ind].sum()
    print(f"The best possible score is: {matchingAlgRes}")
    initialScores(ours, theirs, df)
    print("Who did we pick first? \n")
    ourFirstPick = getFirstPlayer(ours)
    print("Who did they pick first? \n")
    theirFirstPick = getFirstPlayer(theirs)
    printAverageScores(ourFirstPick, theirFirstPick, df)
    return

if __name__ == "__main__":
    main()
