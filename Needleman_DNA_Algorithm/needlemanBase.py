
def result_organizer(align1, align2):
    result = []
    score = 0
    space = ''
    # Reverse Sequences
    align1 = align1[::-1]
    align2 = align2[::-1]

    for i in range(0 ,len(align1)):
        # If the two alignments are the same, then output the former alignment
        if align1[i] == align2[i]:
            space = space + align1[i]

            score += match_score(align1[i], align2[i])

        # If they are not identical and none of them has a gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i])
            space += ' '

        # If one of them is a gap, output a space in the gap
        elif align1[i] == '-' or align2[i] == '-':
            space += ' '
            score += gap_penalty

    # Score is in the first position within the result list
    result.append(score)
    # Determine which result to return for aligned result
    if align1.find(" "):
        result.append(align2)
    else:
        result.append(align1)

    # Returns list containing: [score, aligned result]
    return result


match_reward = 1
mismatch_penalty = -1
gap_penalty = -2

# Function for determining penalty or reward
def match_score( a, b):
    if a == b:
        return match_reward
    elif a == '-' or b == '-':
        return gap_penalty
    else:
        return mismatch_penalty

def scoreHelper(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

#  Core of Needleman-Wunsch algorithm
def needle(seq1, seq2):
    m = len(seq1)
    n = len(seq2)

    # Generate DP table and traceback path pointer matrix
    score = scoreHelper((m + 1, n + 1))

    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment
    align1 = ''
    align2 = ''
    # Start from the bottom right cell
    i = m
    j = n

    # Tracing until it reaches the left edge
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[i - 1], seq2[j - 1]):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    # Returns the list obtained from result_organizer, which will be [score, aligned result]
    return result_organizer(align1, align2)
