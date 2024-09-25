from .aminoacid import Aminoacid, L1, L3, Classification 

class FastaEntry:
    def __init__(self, name: str, sequence: str):
        self.name = name
        self.sequence = sequence

class FastaService:
    def complete_sequences_with_dash(self, fasta_entries: list[FastaEntry]) -> tuple[list[FastaEntry], int]:
        if not fasta_entries:
            return [], 0

        # Encontrar o comprimento máximo das sequências
        max_length = max(len(entry.sequence) for entry in fasta_entries)
        fasta_entries_aux = []

        for entry in fasta_entries:
            diff_length = max_length - len(entry.sequence)
            if diff_length > 0:
                entry.sequence += '.' * diff_length  # Preencher com pontos (ou dashes) para igualar o tamanho
            fasta_entries_aux.append(entry)

        return fasta_entries_aux, max_length

class ListProcessingService:
    def __init__(self):
        pass

    def find_most_frequent_element(self, lst):
        """
        Encontra o elemento mais frequente em uma lista.
        
        Parâmetros:
        - lst: Uma lista de elementos (pode ser de qualquer tipo).
        
        Retorna:
        - O elemento mais frequente ou None se a lista estiver vazia.
        """
        if not lst:  # Verifica se a lista está vazia
            return None

        frequency_map = {}
        most_frequent_element = lst[0]
        max_count = 0

        for element in lst:
            key = str(element)  # Converte o elemento para string para usar como chave no dicionário
            if key in frequency_map:
                frequency_map[key] += 1
            else:
                frequency_map[key] = 1

            # Atualiza o elemento mais frequente se necessário
            if frequency_map[key] > max_count:
                max_count = frequency_map[key]
                most_frequent_element = element

        return most_frequent_element

class PROSITEProcessingService:
    def __init__(self, fasta_service, list_processing_service):
       
        self.fasta_service = fasta_service
        self.list_processing_service = list_processing_service
        self.fasta_entries_with_dashes = []
        self.max_length = 0

    def x_gap_comparate(self, a: str, b: str) -> bool:
        """
        Compara duas strings e retorna True se 'x' for equivalente a 'x0'.
        """
        if (a == 'x' and b == 'x0') or (b == 'x' and a == 'x0'):
            return True
        return False

    def parse_fasta(self, content: str) -> list:
        """
        Analisa o conteúdo de uma sequência FASTA e retorna uma lista de dicionários contendo
        o nome e a sequência de cada entrada.
        """
        lines = content.split('\n')
        entries = []
        name = ''
        sequence = ''

        for line in lines:
            if line.startswith('>'):
                if name and sequence:
                    # Adiciona a entrada anterior antes de iniciar uma nova
                    entries.append({'name': name, 'sequence': sequence})
                name = line[1:].strip()  # Remove o ">" e os espaços em branco
                sequence = ''  # Reseta a sequência para a próxima entrada
            else:
                sequence += line.strip()

        # Adiciona a última entrada, se houver
        if name and sequence:
            entries.append({'name': name, 'sequence': sequence})

        return entries
    
    def scattered_conservation_pattern(self, fasta_entries, score_model_conservation):
        """
        Gera um padrão de conservação dispersa a partir de uma lista de entradas FASTA.
        
        Parâmetros:
        - fasta_entries: Uma lista de objetos FastaEntry.
        - score_model_conservation: O modelo de conservação a ser usado (ex: 'BLOSUM62').
        
        Retorna:
        - Uma lista de strings representando o padrão de conservação dispersa.
        """
        scattered_conservation_pattern = []
        # Completa as sequências com '-' (gaps)
        self.fasta_entries_with_dashes, self.max_length = self.fasta_service.complete_sequences_with_dash(fasta_entries)

        for i in range(self.max_length):
            # Extrai os caracteres de todas as sequências na posição i
            characters_at_position_i = [entry.sequence[i] for entry in self.fasta_entries_with_dashes]
            
            # Encontra o caractere mais frequente na posição atual
            currently_char = self.list_processing_service.find_most_frequent_element(characters_at_position_i)

            # Verifica se todos os caracteres na posição são iguais (conservação completa)
            full_conservation = all(aminoacid == characters_at_position_i[0] for aminoacid in characters_at_position_i)

            # Verifica a conservação com base no modelo de pontuação fornecido
            if score_model_conservation == 'BLOSUM62':
                conservation = all(
                    Aminoacid.get_aminoacid(amino).equals_substitution_matrix_blossum62_1(
                        Aminoacid.get_aminoacid(currently_char)
                    ) for amino in characters_at_position_i
                )
            else:
                conservation = all(
                    Aminoacid.get_aminoacid(amino).equals_classification(
                        Aminoacid.get_aminoacid(currently_char)
                    ) for amino in characters_at_position_i
                )

            if full_conservation:
                if currently_char in ['-', '.']:
                    scattered_conservation_pattern.append('-')
                else:
                    scattered_conservation_pattern.append(currently_char or '0')
            else:
                if conservation:
                    unique_characters = ''.join(set(characters_at_position_i))
                    scattered_conservation_pattern.append(f'[{unique_characters}]')
                else:
                    if '-' in characters_at_position_i:
                        scattered_conservation_pattern.append('x0')
                    else:
                        scattered_conservation_pattern.append('x')

        return scattered_conservation_pattern
    
    def x_min_max_on_the_gaps(self, a: int, b: int) -> (int, int):
        """
        Calcula o valor mínimo e máximo de caracteres diferentes de gaps ('-')
        em um intervalo específico de sequência.

        Parâmetros:
        - a: Início do intervalo (índice).
        - b: Fim do intervalo (índice).

        Retorna:
        - Uma tupla com o valor mínimo e máximo de caracteres não gaps no intervalo.
        """
        result = (0, 0)
        length = b - a + 1
        min_x = float('inf')  # Inicializa o mínimo com infinito
        max_x = 0

        for fasta_entry in self.fasta_entries_with_dashes:
            # Extrai o substring da sequência no intervalo [a, b)
            sub_sequence = fasta_entry.sequence[a:b]

            # Calcula o número de caracteres que não são gaps ('-')
            number_of_x = length - sub_sequence.count('-')

            # Atualiza o valor mínimo
            if number_of_x < min_x:
                min_x = number_of_x - 1

            # Atualiza o valor máximo
            if number_of_x > max_x:
                max_x = number_of_x

        # Se min_x for infinito, significa que todas as posições são gaps, ajusta para zero
        if min_x == float('inf'):
            min_x = 0

        result = (min_x, max_x)
        return result
    
    def group_repeated_strings(self, scattered_conservation_pattern: list[str]) -> list[tuple[list[str], tuple[int, int]]]:
        result = []
        result_aux = []
        range_positions = []
        temp = []

        for i in range(len(scattered_conservation_pattern)):
            current = scattered_conservation_pattern[i]
            
            if self.x_gap_compare(temp[-1] if temp else None, current):
                temp.append('x0')
            elif not temp or temp[-1] == current:
                temp.append(current)
            else:
                result_aux.append(temp.copy())
                temp = [current]

                if len(result_aux) == 1:
                    range_positions.append((0, i))
                    result.append((result_aux[0], range_positions[0]))
                else:
                    range_positions.append((range_positions[-1][1] + 1, i))
                    result.append((result_aux[-1], range_positions[-1]))

        if temp:
            result.append((temp.copy(), (range_positions[-1][1], len(scattered_conservation_pattern) - 1)))

        return result
    
    def count_group_repeated_strings(
        self, 
        fasta_entries: list[FastaEntry], 
        group_repeated_strings: list[tuple[list[str], tuple[int, int]]]
    ) -> list[tuple[str, int, tuple[int, int]]]:
        
        # result será uma lista que armazena os resultados finais
        result = []

        for repeats in group_repeated_strings:
            # Para cada grupo de strings repetidas, adicionamos uma tupla ao resultado
            result.append((repeats[0][0], len(repeats[0]), repeats[1]))

        # Se o primeiro elemento for '-' ou 'x', remove-o
        if result and (result[0][0] == '-' or result[0][0] == 'x'):
            result.pop(0)

        # Se o último elemento for '-' ou 'x', remove-o
        if result and (result[-1][0] == '-' or result[-1][0] == 'x'):
            result.pop()

        return result
    
    def x_threshold_divider(
        self, 
        count_group_repeated_strings: list[tuple[str, int, tuple[int, int]]], 
        xthreshold: int = 20
    ) -> list[list[tuple[str, int, tuple[int, int]]]]:
        
        result = []  # Lista para armazenar os resultados
        aux = []  # Lista auxiliar para agrupar strings

        for group in count_group_repeated_strings:
            # Verifica a condição para adicionar ao grupo auxiliar
            if ((group[0] == 'x' or group[0] == 'x0') and group[1] < xthreshold) or (group[0] != 'x' and group[0] != 'x0'):
                aux.append(group)
            elif aux:  # Se a lista auxiliar não estiver vazia
                result.append(aux)
                aux = []  # Limpa a lista auxiliar

        if aux:  # Se ainda houver elementos na lista auxiliar, adicione ao resultado
            result.append(aux)

        return result
    
    def fitx_threshold_divider(
        self, 
        PROSITE_motifs_pattern_x_threshold_divided: list[list[tuple[str, int, tuple[int, int]]]]
    ) -> list[list[tuple[str, int, tuple[int, int]]]]:
        
        result = []  # Lista para armazenar os resultados
        
        for motif in PROSITE_motifs_pattern_x_threshold_divided:
            if motif[0][0] == '-':
                motif.pop(0)  # Remove o primeiro elemento se for '-'
                
            # Se o comprimento do motif não for zero, adiciona ao resultado
            if motif:
                result.append(motif)

        print(result)  # Usar print para depuração (equivalente ao console.log em JS)
        return result
    
    def format_prosite_motifs_pattern(
        self, 
        PROSITE_motifs_pattern_x_threshold_divided: list[list[tuple[str, int, tuple[int, int]]]]
    ) -> list[list[str]]:
        
        result = []  # Lista para armazenar os resultados
        aux = []  # Lista auxiliar para armazenar os padrões formatados

        for pattern in PROSITE_motifs_pattern_x_threshold_divided:
            for conservation in pattern:
                if conservation[1] > 1:
                    if conservation[0] == 'x0':
                        auxx = self.x_min_max_on_the_gaps(conservation[2][0], conservation[2][1])
                        aux.append(f'x({auxx[0]},{auxx[1]})')  # Formata 'x(xmin,xmax)'
                    else:
                        aux.append(f'{conservation[0]}({conservation[1]})')  # Formata 'element(count)'
                else:
                    aux.append(conservation[0])  # Apenas adiciona o elemento

            result.append(aux)  # Adiciona o padrão formatado à lista de resultados
            aux = []  # Reseta a lista auxiliar para o próximo padrão

        return result  # Retorna a lista de resultados
    
    def process_fasta(
        self, 
        fasta_content: str, 
        score_model_conservation: str, 
        xthreshold: int | None
    ) -> list[list[str]]:
        
        # Parse the FASTA content into entries
        fasta_entries = self.parse_fasta(fasta_content)

        # Process the parsed entries through the defined methods
        return self.format_prosite_motifs_pattern(
            self.fit_x_treshold_divider(
                self.x_treshold_divider(
                    self.count_group_repeated_strings(
                        fasta_entries,
                        self.group_repeated_strings(
                            self.scattered_conservation_pattern(
                                fasta_entries,
                                score_model_conservation
                            )
                        )
                    ),
                    xthreshold
                )
            )
        )